// Shiny output binding that draws a Plotly figure sent from the server as plain JSON
// ({data, layout, config}) with plotly.js. Used instead of shinywidgets, which
// re-sent the whole plotly.js library (~5 MB) with every render.
(function () {
  class PlotlyOutputBinding extends Shiny.OutputBinding {
    find(scope) {
      return $(scope).find(".plotly-output");
    }

    renderValue(el, fig) {
      if (!fig) {
        this._clear(el);
        return;
      }
      let plot = el.querySelector(":scope > .plotly-target");
      if (!plot) {
        el.replaceChildren();
        plot = document.createElement("div");
        plot.className = "plotly-target";
        el.appendChild(plot);
        this._watchSize(plot);
      }
      // Without a fixed height, fill the output container (as shinywidgets did)
      plot.classList.toggle("html-fill-item", fig.layout.height == null);
      // Draw from scratch each time (not Plotly.react), as shinywidgets did, so the
      // box plots' point jitter is laid out the same way on every render
      Plotly.purge(plot);
      Plotly.newPlot(plot, fig.data, fig.layout, fig.config);
    }

    renderError(el, err) {
      this._clear(el);
      super.renderError(el, err);
    }

    // A plot drawn while its tab is hidden can't measure its width, so plotly falls
    // back to 700px, and fixed-height plots never resize after that. Redraw at the
    // container's size whenever it changes, including when the tab is shown.
    _watchSize(plot) {
      let pending = false;
      plot._sizeObserver = new ResizeObserver(() => {
        if (pending) return;
        pending = true;
        requestAnimationFrame(() => {
          pending = false;
          if (plot._fullLayout && plot.offsetWidth > 0) Plotly.Plots.resize(plot);
        });
      });
      plot._sizeObserver.observe(plot);
    }

    _clear(el) {
      const plot = el.querySelector(":scope > .plotly-target");
      if (plot) {
        plot._sizeObserver?.disconnect();
        Plotly.purge(plot);
      }
      el.replaceChildren();
    }
  }

  Shiny.outputBindings.register(new PlotlyOutputBinding(), "aging.plotlyOutput");
})();
