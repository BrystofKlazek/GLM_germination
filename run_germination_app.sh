#!/usr/bin/env bash
# ──────────────────────────────────────────────────────────────
#  GLM Germination Analysis — Shiny App Launcher (macOS / Linux)
# ──────────────────────────────────────────────────────────────
set -euo pipefail

echo ""
echo "  ============================================================"
echo "    GLM Germination Analysis — Shiny App Launcher"
echo "  ============================================================"
echo ""

# ── Locate Rscript ──────────────────────────────────────────────
RSCRIPT=""

if command -v Rscript &>/dev/null; then
    RSCRIPT="$(command -v Rscript)"
fi

if [ -z "$RSCRIPT" ]; then
    for candidate in \
        /usr/local/bin/Rscript \
        /opt/homebrew/bin/Rscript \
        /Library/Frameworks/R.framework/Resources/bin/Rscript \
        /usr/bin/Rscript \
        /usr/lib/R/bin/Rscript; do
        if [ -x "$candidate" ] 2>/dev/null; then
            RSCRIPT="$candidate"
            break
        fi
    done
fi

if [ -z "$RSCRIPT" ]; then
    echo "  [ERROR] Rscript not found!"
    echo "    Install R: sudo pacman -S r gcc-fortran  (Arch/CachyOS)"
    echo "               brew install r                (macOS)"
    echo "               sudo apt install r-base       (Debian/Ubuntu)"
    exit 1
fi

echo "  [OK] Found R: $RSCRIPT"
R_VERSION=$("$RSCRIPT" --version 2>&1 | head -1)
echo "       $R_VERSION"
echo ""

# ── Locate app.R ────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
APP_FILE="$SCRIPT_DIR/app.R"

if [ ! -f "$APP_FILE" ]; then
    echo "  [ERROR] app.R not found in $SCRIPT_DIR"
    exit 1
fi

echo "  [OK] App found: $APP_FILE"
echo ""

# ── Ensure BROWSER is set ──────────────────────────────────────
if [ -z "${BROWSER:-}" ]; then
    # Try to auto-detect a browser
    for b in xdg-open firefox chromium google-chrome-stable brave; do
        if command -v "$b" &>/dev/null; then
            export BROWSER="$b"
            break
        fi
    done
fi

if [ -z "${BROWSER:-}" ]; then
    echo "  [NOTE] No browser detected. The app will start but you need to"
    echo "         open http://127.0.0.1:PORT manually in your browser."
    echo ""
fi

# ── Check packages ─────────────────────────────────────────────
echo "  Checking required R packages..."
echo ""

"$RSCRIPT" --no-save --no-restore -e '
pkgs <- c("shiny","bslib","readxl","readr","emmeans","ggplot2",
          "dplyr","tidyr","broom","DT","multcomp","car","scales")
missing <- setdiff(pkgs, rownames(installed.packages()))
if (length(missing)) {
  cat("  Installing:", paste(missing, collapse = ", "), "\n")
  install.packages(missing, repos = "https://cloud.r-project.org", quiet = TRUE)
}
cat("  All packages OK.\n")
' || {
    echo ""
    echo "  [WARNING] Package check had issues. The app will try on start."
    echo ""
}

# ── Launch ────────────────────────────────────────────────────
echo ""
echo "  ============================================================"
echo "    Launching Shiny app..."
if [ -n "${BROWSER:-}" ]; then
    echo "    Browser: $BROWSER"
else
    echo "    Open http://127.0.0.1:PORT in your browser manually."
fi
echo "    Press Ctrl+C to stop the app."
echo "  ============================================================"
echo ""

"$RSCRIPT" --no-save --no-restore -e "
browser_cmd <- Sys.getenv('BROWSER', '')
if (nzchar(browser_cmd)) {
  options(browser = browser_cmd)
  shiny::runApp('$APP_FILE', launch.browser = TRUE)
} else {
  cat('  Open this URL in your browser: ')
  shiny::runApp('$APP_FILE', launch.browser = function(url) cat(url, '\n'))
}
"

echo ""
echo "  App stopped."
