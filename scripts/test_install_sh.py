#!/usr/bin/env python3
"""Tests for the checksum verification of install.sh.

The installer is sourced without its final ``main "$@"`` line and its functions
are called in a subprocess. ``curl`` is replaced by a fake that serves files from
a local directory, so nothing touches the network or a system path.
"""

import hashlib
import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
INSTALL_SH = ROOT.joinpath("install.sh")

ASSET = "zsasa-0.9.1-linux-x86_64"
BINARY = b"#!/bin/sh\necho zsasa 0.9.1\n"
BINARY_SHA = hashlib.sha256(BINARY).hexdigest()

# Tools install.sh needs besides curl and a sha256 tool.
BASE_TOOLS = (
    "awk",
    "basename",
    "cat",
    "cp",
    "install",
    "mkdir",
    "mktemp",
    "rm",
    "uname",
)
FAKE_CURL = """\
#!/bin/sh
# Imitates `curl ... -o DEST URL` with files from $FAKE_RELEASE_DIR; a missing file is an HTTP error (exit 22).
dest=""
while [ $# -gt 1 ]; do
    if [ "$1" = "-o" ]; then dest="$2"; shift; fi
    shift
done
src="$FAKE_RELEASE_DIR/$(basename "$1")"
[ -f "$src" ] || exit 22
cp "$src" "$dest"
"""


class InstallShTests(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="zsasa-install-test-"))
        self.addCleanup(lambda: shutil.rmtree(self.tmp, ignore_errors=True))
        self.release = self.tmp.joinpath("release")
        self.release.mkdir()
        self.install_dir = self.tmp.joinpath("install")
        self.shell = os.environ.get("INSTALL_SH_TEST_SHELL") or shutil.which(
            "sh"
        )  # e.g. dash, to check POSIX conformance
        lines = INSTALL_SH.read_text().rstrip("\n").split("\n")
        self.assertEqual(
            lines[-1], 'main "$@"', "install.sh must end with the call to main"
        )
        self.lib = self.tmp.joinpath("install_lib.sh")
        self.lib.write_text("\n".join(lines[:-1]) + "\n")

    def make_bin(self, sha_tools: tuple[str, ...]) -> Path:
        """Return a PATH directory with the fake curl, the basic tools and the given sha256 tools only."""
        bin_dir = self.tmp.joinpath("bin-" + "-".join(sha_tools or ("none",)))
        if bin_dir.exists():
            return bin_dir
        bin_dir.mkdir()
        for tool in (*BASE_TOOLS, *sha_tools):
            real = shutil.which(tool)
            self.assertIsNotNone(real, f"{tool} not found on this machine")
            bin_dir.joinpath(tool).symlink_to(real)
        curl = bin_dir.joinpath("curl")
        curl.write_text(FAKE_CURL)
        curl.chmod(0o755)
        return bin_dir

    def call(
        self,
        script: str,
        *,
        sha_tools: tuple[str, ...] = ("sha256sum",),
        env: dict[str, str] | None = None,
    ):
        if "sha256sum" in sha_tools and shutil.which("sha256sum") is None:
            sha_tools = ("shasum",)  # macOS ships shasum only
        environment = {
            "PATH": str(self.make_bin(sha_tools)),
            "HOME": str(self.tmp),
            "FAKE_RELEASE_DIR": str(self.release),
            **(env or {}),
        }
        return subprocess.run(
            [self.shell, "-c", f'. "{self.lib}"\n{script}'],
            env=environment,
            text=True,
            capture_output=True,
            check=False,
        )

    def publish(self, *, sums: str | None, binary: bytes = BINARY):
        """Put the asset (and, when given, the SHA256SUMS text) into the fake release."""
        self.release.joinpath(ASSET).write_bytes(binary)
        sums_path = self.release.joinpath("SHA256SUMS")
        sums_path.unlink(missing_ok=True)
        if sums is not None:
            sums_path.write_text(sums)

    def installed(self) -> bool:
        return self.install_dir.joinpath("zsasa").exists()

    def download(self, **kwargs):
        return self.call(
            f'make_tmpdir\ndownload_binary 0.9.1 linux-x86_64 "{self.install_dir}"',
            **kwargs,
        )

    # -- verify_checksum -------------------------------------------------
    def verify(self, sums: str, **kwargs):
        binary = self.tmp.joinpath("zsasa")
        binary.write_bytes(BINARY)
        sums_file = self.tmp.joinpath("SHA256SUMS")
        sums_file.write_text(sums)
        return self.call(f'verify_checksum "{binary}" {ASSET} "{sums_file}"', **kwargs)

    def test_verify_accepts_the_matching_checksum(self):
        for sums in (
            f"{BINARY_SHA}  {ASSET}\n",
            f"{'0' * 64}  zsasa-0.9.1-macos-aarch64\n{BINARY_SHA}  {ASSET}\n",
            f"{BINARY_SHA} *{ASSET}\n",  # sha256sum --binary
            f"{BINARY_SHA.upper()}  {ASSET}\n",
        ):
            result = self.verify(sums)
            self.assertEqual(result.returncode, 0, (sums, result.stderr))
            self.assertIn("Checksum verified.", result.stdout)

    def test_verify_works_with_shasum_only(self):
        if shutil.which("shasum") is None:
            self.skipTest("shasum not installed")
        result = self.verify(f"{BINARY_SHA}  {ASSET}\n", sha_tools=("shasum",))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("Checksum verified.", result.stdout)

    def test_verify_rejects_a_mismatch(self):
        result = self.verify(f"{'1' * 64}  {ASSET}\n")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Checksum mismatch", result.stderr)
        self.assertNotIn("Checksum verified.", result.stdout)

    def test_verify_requires_an_entry_for_exactly_this_asset(self):
        for sums in (
            "",
            "\n",
            f"{BINARY_SHA}  zsasa-0.9.1-macos-aarch64\n",
            f"{BINARY_SHA}  {ASSET}.exe\n",  # a longer name that merely contains the asset name
            f"{BINARY_SHA}  other-{ASSET}\n",
        ):
            result = self.verify(sums)
            self.assertNotEqual(result.returncode, 0, repr(sums))
            self.assertIn("no entry", result.stderr)
            self.assertNotIn("Checksum verified.", result.stdout)

    def test_verify_rejects_malformed_checksums(self):
        for digest in ("abc", "g" * 64, BINARY_SHA + "0"):
            result = self.verify(f"{digest}  {ASSET}\n")
            self.assertNotEqual(result.returncode, 0, digest)
            self.assertIn("malformed", result.stderr)

    def test_verify_fails_without_a_sha256_tool(self):
        result = self.verify(f"{BINARY_SHA}  {ASSET}\n", sha_tools=())
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Neither sha256sum nor shasum", result.stderr)
        self.assertNotIn("Checksum verified.", result.stdout)

    # -- download_binary -------------------------------------------------
    def test_download_installs_after_a_successful_verification(self):
        self.publish(sums=f"{BINARY_SHA}  {ASSET}\n")
        result = self.download()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("Checksum verified.", result.stdout)
        self.assertEqual(self.install_dir.joinpath("zsasa").read_bytes(), BINARY)

    def test_download_stops_when_sha256sums_cannot_be_downloaded(self):
        self.publish(sums=None)
        result = self.download()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Could not download", result.stderr)
        self.assertIn("SKIP_CHECKSUM=1", result.stderr)
        self.assertNotIn("Checksum verified.", result.stdout)
        self.assertFalse(self.installed())

    def test_download_stops_when_sha256sums_has_no_entry(self):
        self.publish(sums=f"{BINARY_SHA}  zsasa-0.9.1-macos-aarch64\n")
        result = self.download()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("no entry", result.stderr)
        self.assertFalse(self.installed())

    def test_download_stops_on_a_checksum_mismatch(self):
        self.publish(sums=f"{BINARY_SHA}  {ASSET}\n", binary=BINARY + b"tampered")
        result = self.download()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Checksum mismatch", result.stderr)
        self.assertFalse(self.installed())

    def test_download_stops_without_a_sha256_tool(self):
        self.publish(sums=f"{BINARY_SHA}  {ASSET}\n")
        result = self.download(sha_tools=())
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Neither sha256sum nor shasum", result.stderr)
        self.assertNotIn("Checksum verified.", result.stdout)
        self.assertFalse(self.installed())

    def test_skip_checksum_installs_without_verifying(self):
        self.publish(sums=None)
        result = self.download(env={"SKIP_CHECKSUM": "1"})
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("SKIP_CHECKSUM=1", result.stderr)
        self.assertNotIn("Checksum verified.", result.stdout)
        self.assertEqual(self.install_dir.joinpath("zsasa").read_bytes(), BINARY)

    def test_skip_checksum_needs_the_value_1(self):
        self.publish(sums=None)
        for value in ("0", "", "yes", "true"):
            result = self.download(env={"SKIP_CHECKSUM": value})
            self.assertNotEqual(result.returncode, 0, value)
            self.assertFalse(self.installed(), value)

    # -- the script itself -----------------------------------------------
    def test_script_is_valid_posix_sh_and_documents_the_opt_out(self):
        subprocess.run([self.shell, "-n", str(INSTALL_SH)], check=True)
        usage = INSTALL_SH.read_text().split("set -eu")[0]
        self.assertIn("SKIP_CHECKSUM", usage)

    def test_installer_has_no_other_way_to_skip_verification(self):
        text = INSTALL_SH.read_text()
        self.assertNotIn("Skipping", text)
        self.assertEqual(text.count("SKIP_CHECKSUM:-"), 1)


if __name__ == "__main__":
    unittest.main()
