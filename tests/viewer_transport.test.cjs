const { test } = require('node:test');
const assert = require('node:assert/strict');
require('../src/paratan/viewer/static/transport.js');
const { deduplicateArrays } = globalThis.paratan_transport;

test('concurrent geometry reads share one request; later reads use the VTK cache', async () => {
  let calls = 0;
  let complete;
  const array = new Float32Array([1, 2, 3]);
  const context = { getArray: () => { calls++; return new Promise(resolve => { complete = resolve; }); } };
  deduplicateArrays(context);
  deduplicateArrays(context);
  const reads = Array.from({ length: 50 }, () => context.getArray('mesh', 'Float32Array', {}));
  await Promise.resolve();
  assert.equal(calls, 1);
  complete(array);
  assert.ok((await Promise.all(reads)).every(value => value === array));
  const next = context.getArray('mesh', 'Float32Array', {});
  await Promise.resolve();
  assert.equal(calls, 2);
  complete(array);
  await next;
});

test('failed reads can be retried and different types are independent', async () => {
  let calls = 0;
  const context = { getArray: async () => { if (++calls === 1) throw new Error('network'); return calls; } };
  deduplicateArrays(context);
  await assert.rejects(context.getArray('mesh', 'Float32Array'), /network/);
  await context.getArray('mesh', 'Float32Array');
  await Promise.all([context.getArray('mesh', 'Float32Array'), context.getArray('mesh', 'Uint32Array')]);
  assert.equal(calls, 4);
});
