/* Share concurrent VTK array reads without changing VTK's retained array cache. */
(function (root) {
  function deduplicateArrays(context) {
    if (context.paratanSharedReads) return;
    const read = context.getArray.bind(context);
    const pending = new Map();
    context.getArray = function (hash, type, ...rest) {
      const key = `${hash}:${type}`;
      if (!pending.has(key)) {
        const request = Promise.resolve().then(() => read(hash, type, ...rest));
        pending.set(key, request.finally(() => pending.delete(key)));
      }
      return pending.get(key);
    };
    context.paratanSharedReads = true;
  }

  root.paratan_transport = {
    deduplicateArrays,
    install(app) {
      app.component('paratan-scene-transport', {
        setup() {
          const view = root.Vue.inject('view', null);
          if (view && view.ctx) deduplicateArrays(view.ctx);
          return () => null;
        },
      });
    },
  };
})(typeof window === 'undefined' ? globalThis : window);
