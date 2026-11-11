(() => {
  function markExternalLinks() {
    for (const link of document.querySelectorAll("a[href]")) {
      const url = new URL(link.href, window.location.href);
      const isWebLink = url.protocol === "http:" || url.protocol === "https:";

      if (isWebLink && url.origin !== window.location.origin) {
        link.target = "_blank";
        link.rel = "noopener noreferrer";
      }
    }
  }

  document.addEventListener("DOMContentLoaded", markExternalLinks);

  if (typeof document$ !== "undefined") {
    document$.subscribe(markExternalLinks);
  }
})();
