/* Drag-and-drop for #paratan-upload-dropzone — keeps Vue templates free of drop handlers. */
(function () {
  function state() {
    return window.trame && window.trame.state;
  }

  function set(key, value) {
    const s = state();
    if (s && typeof s.set === "function") s.set(key, value);
  }

  function get(key) {
    const s = state();
    return s && typeof s.get === "function" ? s.get(key) : undefined;
  }

  function uploadFile(file) {
    if (!file || get("upload_busy")) return;
    const name = file.name || "input.yaml";
    const lower = name.toLowerCase();
    if (!(lower.endsWith(".yaml") || lower.endsWith(".yml"))) {
      set("upload_phase", "error");
      set("upload_notice", true);
      set("load_status", name + " is not a .yaml / .yml file.");
      return;
    }
    set("upload_drag_over", false);
    set("upload_busy", true);
    set("upload_phase", "loading");
    set("upload_notice", false);
    set("load_status", "Loading " + name + "…");
    if (file.size > 2000000) {
      set("upload_busy", false);
      set("upload_phase", "error");
      set("upload_notice", true);
      set("load_status", name + " is larger than 2 MB.");
      return;
    }
    const aborter = new AbortController();
    const timer = setTimeout(() => aborter.abort(), 30000);
    fetch("/paratan-upload?name=" + encodeURIComponent(name), {
      method: "POST",
      body: file,
      signal: aborter.signal,
    })
      .then((r) => r.json())
      .then((result) => {
        set("upload_phase", result.phase);
        set("load_status", result.message);
        set("upload_notice", true);
      })
      .catch((error) => {
        set("upload_phase", "error");
        set("upload_notice", true);
        set(
          "load_status",
          "Could not upload " + name + ": " + error.message + ". Please try again."
        );
      })
      .finally(() => {
        clearTimeout(timer);
        set("upload_busy", false);
      });
  }

  function bind(el) {
    if (!el || el.__paratanDropBound) return;
    el.__paratanDropBound = true;
    let depth = 0;
    function highlight(on) {
      el.classList.toggle("paratan-dropzone--active", on);
      set("upload_drag_over", on);
    }
    el.addEventListener("dragenter", (e) => {
      e.preventDefault();
      if (get("upload_busy")) return;
      depth += 1;
      highlight(true);
    });
    el.addEventListener("dragover", (e) => {
      e.preventDefault();
      if (get("upload_busy")) return;
      highlight(true);
    });
    el.addEventListener("dragleave", (e) => {
      e.preventDefault();
      depth = Math.max(0, depth - 1);
      if (!depth) highlight(false);
    });
    el.addEventListener("drop", (e) => {
      e.preventDefault();
      depth = 0;
      highlight(false);
      const file = e.dataTransfer && e.dataTransfer.files && e.dataTransfer.files[0];
      uploadFile(file);
    });
  }

  function tryBind() {
    bind(document.getElementById("paratan-upload-dropzone"));
  }

  // Trame mounts asynchronously; retry briefly until the dropzone exists.
  let tries = 0;
  const timer = setInterval(() => {
    tryBind();
    tries += 1;
    if (document.getElementById("paratan-upload-dropzone") || tries > 40) {
      clearInterval(timer);
    }
  }, 250);
  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", tryBind);
  } else {
    tryBind();
  }
})();
