# Phase 14: WebAssembly Support

## Goal

WebAssembly (WASM) にコンパイルし、ブラウザやNode.jsから利用可能にする。

## 優先度: 極低

## 概要

ZigはWASMターゲットをネイティブサポートしているため、比較的容易にWASMビルドが可能。

## Tasks

### Phase 14.1: WASM Build

- [ ] wasm32-freestanding ターゲットビルド
- [ ] JavaScript glue code生成
- [ ] メモリ管理（WASMヒープ）

**ビルド設定:**
```zig
// build.zig
const wasm = b.addSharedLibrary(.{
    .name = "freesasa_zig",
    .root_source_file = .{ .path = "src/wasm.zig" },
    .target = .{ .cpu_arch = .wasm32, .os_tag = .freestanding },
    .optimize = .ReleaseSmall,
});
```

### Phase 14.2: JavaScript API

- [ ] TypeScript型定義
- [ ] 非同期API（Web Worker対応）
- [ ] ストリーミング入力

**JavaScript API:**
```typescript
import init, { calculateSasa } from 'freesasa-zig';

await init();

const result = calculateSasa({
  x: Float64Array,
  y: Float64Array,
  z: Float64Array,
  r: Float64Array,
  nPoints: 100,
  probeRadius: 1.4,
});

console.log(result.totalArea);
console.log(result.atomAreas);
```

### Phase 14.3: npm Package

- [ ] npm パッケージング
- [ ] ESM/CommonJS両対応
- [ ] ドキュメント

## Files

| File | Action |
|------|--------|
| `src/wasm.zig` | CREATE |
| `wasm/index.js` | CREATE |
| `wasm/index.d.ts` | CREATE |
| `wasm/package.json` | CREATE |
| `build.zig` | MODIFY |

## Limitations

- マルチスレッド: Web Workers経由で可能だが複雑
- SIMD: WASM SIMDサポートが必要（ブラウザ依存）
- メモリ: 大規模構造では制限あり

## Success Criteria

- [ ] ブラウザで動作
- [ ] Node.jsで動作
- [ ] npmでインストール可能

---
- [ ] **DONE** - Phase 14 complete
