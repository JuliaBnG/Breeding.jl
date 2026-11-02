// JuliaBnG Sepia and Black themes for Documenter.
//
// Documenter's theme list is fixed, so these themes come in as extra assets:
// <link data-theme-name="juliabng-..."> stylesheets placed after themeswap.js,
// then this script. Here we re-run Documenter's swap now that those sheets
// exist, and add the two themes to the Settings theme picker. Documenter's own
// picker handler then stores and applies them like any built-in theme
// (localStorage "documenter-theme").
(function () {
  var THEMES = [["juliabng-sepia", "sepia"], ["juliabng-black", "black"]];

  if (typeof set_theme_from_local_storage === "function") {
    set_theme_from_local_storage();
  }

  function addOptions() {
    var sel = document.getElementById("documenter-themepicker");
    if (!sel) return;
    var after = sel.querySelector('option[value="documenter-dark"]');
    THEMES.forEach(function (t) {
      if (sel.querySelector('option[value="' + t[0] + '"]')) return;
      var o = document.createElement("option");
      o.value = t[0];
      o.textContent = t[1];
      if (after) { after.after(o); after = o; } else { sel.appendChild(o); }
    });
    // Documenter syncs the picker on its own ready handler, which may have run
    // before these options existed.
    var cur = null;
    try { cur = window.localStorage.getItem("documenter-theme"); } catch (e) {}
    if (cur && sel.querySelector('option[value="' + cur + '"]')) sel.value = cur;
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", addOptions);
  } else {
    addOptions();
  }
})();
