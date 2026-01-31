# Security Policy

## Supported Versions

| Version | Supported          |
| ------- | ------------------ |
| 0.1.x   | :white_check_mark: |
| < 0.1.0 | :x:                |

## Reporting a Vulnerability

If you discover a security vulnerability in freesasa-zig, please report it privately:

1. **GitHub Security Advisories** (preferred): Use [GitHub's private vulnerability reporting](https://github.com/N283T/freesasa-zig/security/advisories/new)
2. **Email**: Contact the maintainer directly

Please include:
- Description of the vulnerability
- Steps to reproduce
- Potential impact
- Suggested fix (if any)

We aim to respond within 48 hours and will work with you to understand and address the issue.

## Security Considerations

freesasa-zig processes molecular structure data. While the library itself does not handle sensitive data, users should be aware of:

- **Input validation**: The library validates input data (coordinates, radii) but malformed files could potentially cause issues
- **Memory safety**: Written in Zig with explicit memory management; no garbage collection
- **No network access**: The library does not make network connections
- **No file system writes**: Unless explicitly requested via CLI output arguments
