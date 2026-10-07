(() => {
  const root = document.documentElement;

  // Theme toggle. The initial theme is applied by an inline script in <head>.
  document.querySelector("[data-theme-toggle]")?.addEventListener("click", () => {
    const current =
      root.dataset.theme || "dark";
    const next = current === "dark" ? "light" : "dark";
    root.dataset.theme = next;
    try {
      localStorage.setItem("zsasa-theme", next);
    } catch {}
    window.dispatchEvent(new Event("themechange"));
  });

  // Tabs.
  for (const el of document.querySelectorAll("[data-tabs]")) {
    const tabs = [...el.querySelectorAll(":scope > .tabs__list > [role=tab]")];
    const panels = [...el.querySelectorAll(":scope > .tabs__panel")];
    const select = (i, focus) => {
      tabs.forEach((t, k) => {
        t.setAttribute("aria-selected", String(k === i));
        t.tabIndex = k === i ? 0 : -1;
      });
      panels.forEach((p, k) => (p.hidden = k !== i));
      if (focus) tabs[i].focus();
    };
    tabs.forEach((t, i) => {
      t.addEventListener("click", () => select(i));
      t.addEventListener("keydown", (e) => {
        if (e.key === "ArrowRight") select((i + 1) % tabs.length, true);
        if (e.key === "ArrowLeft") select((i - 1 + tabs.length) % tabs.length, true);
      });
    });
  }

  // Copy buttons.
  for (const btn of document.querySelectorAll("[data-copy]")) {
    btn.addEventListener("click", async () => {
      const pre = btn.closest(".code")?.querySelector("pre");
      if (!pre) return;
      try {
        await navigator.clipboard.writeText(pre.innerText.replace(/\n$/, ""));
        btn.textContent = "Copied";
      } catch {
        btn.textContent = "Failed";
      }
      setTimeout(() => (btn.textContent = "Copy"), 1400);
    });
  }

  // Reveal-on-scroll (bar chart).
  const reveal = document.querySelectorAll("[data-reveal]");
  if (reveal.length) {
    const io = new IntersectionObserver(
      (entries) => {
        for (const e of entries) {
          if (e.isIntersecting) {
            e.target.classList.add("in");
            io.unobserve(e.target);
          }
        }
      },
      { threshold: 0.25 },
    );
    reveal.forEach((el) => io.observe(el));
  }

  // Docs: mobile sidebar.
  const menu = document.querySelector("[data-menu]");
  menu?.addEventListener("click", () => {
    const open = document.body.classList.toggle("menu-open");
    menu.setAttribute("aria-expanded", String(open));
  });

  // Docs: table-of-contents scroll spy.
  const links = [...document.querySelectorAll(".toc a")];
  if (links.length) {
    const heads = links.map((a) => document.getElementById(decodeURIComponent(a.hash.slice(1))));
    const spy = () => {
      let on = 0;
      heads.forEach((h, i) => {
        if (h && h.getBoundingClientRect().top < 120) on = i;
      });
      links.forEach((a, i) => a.classList.toggle("on", i === on));
    };
    document.addEventListener("scroll", spy, { passive: true });
    spy();
  }
})();
