# Naming Convention

## Current

| Item | Name | Convention |
|------|------|------------|
| Repository | `freesasa-zig` | kebab-case (GitHub standard) |
| Binary | `freesasa_zig` | snake_case |
| Files | `shrake_rupley.zig` | snake_case |
| Functions | `calculateSasa` | camelCase |

## Future: Package Release

When releasing as a Zig package, consider renaming to PascalCase:

```
FreeSasaZig
```

Zig packages that export types directly typically use PascalCase.

## Reference

- [Zig naming conventions](https://nathancraddock.com/blog/zig-naming-conventions/)
