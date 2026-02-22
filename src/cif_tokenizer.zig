//! CIF (Crystallographic Information File) Tokenizer.
//!
//! This module provides a simple tokenizer for parsing CIF/mmCIF files.
//! It handles the basic CIF syntax elements needed for atom_site extraction.
//!
//! CIF Syntax Elements:
//! - `data_<name>` - Data block declaration
//! - `loop_` - Loop block declaration
//! - `_category.field` - Tag names
//! - Values: unquoted, single-quoted, double-quoted, or semicolon-delimited
//! - Comments: lines starting with #
//!
//! Reference: https://www.iucr.org/resources/cif/spec/version1.1/cifsyntax

const std = @import("std");

/// Token types for CIF parsing
pub const Token = union(enum) {
    /// Data block declaration: data_<name>
    data_block: []const u8,
    /// Loop declaration: loop_
    loop,
    /// Save frame begin: save_<name> (not commonly used in mmCIF)
    save_begin: []const u8,
    /// Save frame end: save_
    save_end,
    /// Tag name: _category.field
    tag: []const u8,
    /// Value (string, number, or special value like . or ?)
    value: []const u8,
    /// End of file
    eof,
    /// Error token with message
    err: []const u8,

    /// Check if this is a value token
    pub fn isValue(self: Token) bool {
        return self == .value;
    }

    /// Check if this is a tag token
    pub fn isTag(self: Token) bool {
        return self == .tag;
    }

    /// Get the string value if this is a value or tag token
    pub fn getString(self: Token) ?[]const u8 {
        return switch (self) {
            .value => |v| v,
            .tag => |t| t,
            .data_block => |d| d,
            .save_begin => |s| s,
            else => null,
        };
    }
};

/// CIF Tokenizer
pub const Tokenizer = struct {
    source: []const u8,
    pos: usize,
    line: usize,
    col: usize,

    /// Initialize tokenizer with source text
    pub fn init(source: []const u8) Tokenizer {
        return .{
            .source = source,
            .pos = 0,
            .line = 1,
            .col = 1,
        };
    }

    /// Get the next token
    pub fn next(self: *Tokenizer) Token {
        self.skipWhitespaceAndComments();

        if (self.pos >= self.source.len) {
            return .eof;
        }

        const c = self.source[self.pos];

        // Check for semicolon text block (must be at beginning of line)
        if (c == ';' and self.col == 1) {
            return self.readSemicolonText();
        }

        // Check for quoted strings
        if (c == '\'' or c == '"') {
            return self.readQuotedString(c);
        }

        // Check for keywords and tags
        if (c == '_') {
            return self.readTag();
        }

        if (c == 'd' or c == 'D') {
            if (self.matchKeyword("data_")) {
                return self.readDataBlock();
            }
        }

        if (c == 'l' or c == 'L') {
            if (self.matchKeyword("loop_")) {
                self.pos += 5;
                self.col += 5;
                return .loop;
            }
        }

        if (c == 's' or c == 'S') {
            if (self.matchKeyword("save_")) {
                return self.readSaveFrame();
            }
        }

        // Otherwise, read as unquoted value
        return self.readUnquotedValue();
    }

    /// Skip whitespace and comments
    fn skipWhitespaceAndComments(self: *Tokenizer) void {
        while (self.pos < self.source.len) {
            const c = self.source[self.pos];

            if (c == ' ' or c == '\t' or c == '\r') {
                self.pos += 1;
                self.col += 1;
            } else if (c == '\n') {
                self.pos += 1;
                self.line += 1;
                self.col = 1;
            } else if (c == '#') {
                // Skip comment until end of line
                while (self.pos < self.source.len and self.source[self.pos] != '\n') {
                    self.pos += 1;
                }
            } else {
                break;
            }
        }
    }

    /// Check if remaining source matches a keyword (case-insensitive)
    fn matchKeyword(self: *Tokenizer, keyword: []const u8) bool {
        if (self.pos + keyword.len > self.source.len) {
            return false;
        }

        for (keyword, 0..) |kc, i| {
            const sc = self.source[self.pos + i];
            if (std.ascii.toLower(sc) != std.ascii.toLower(kc)) {
                return false;
            }
        }
        return true;
    }

    /// Read a data block name: data_<name>
    fn readDataBlock(self: *Tokenizer) Token {
        self.pos += 5; // skip "data_"
        self.col += 5;

        const start = self.pos;
        while (self.pos < self.source.len) {
            const c = self.source[self.pos];
            if (isWhitespace(c)) break;
            self.pos += 1;
            self.col += 1;
        }

        return .{ .data_block = self.source[start..self.pos] };
    }

    /// Read a save frame: save_<name> or save_ (end)
    fn readSaveFrame(self: *Tokenizer) Token {
        self.pos += 5; // skip "save_"
        self.col += 5;

        const start = self.pos;
        while (self.pos < self.source.len) {
            const c = self.source[self.pos];
            if (isWhitespace(c)) break;
            self.pos += 1;
            self.col += 1;
        }

        if (start == self.pos) {
            return .save_end;
        }
        return .{ .save_begin = self.source[start..self.pos] };
    }

    /// Read a tag: _category.field
    fn readTag(self: *Tokenizer) Token {
        const start = self.pos;
        self.pos += 1; // skip '_'
        self.col += 1;

        while (self.pos < self.source.len) {
            const c = self.source[self.pos];
            if (isWhitespace(c)) break;
            self.pos += 1;
            self.col += 1;
        }

        return .{ .tag = self.source[start..self.pos] };
    }

    /// Read a quoted string (single or double quotes)
    fn readQuotedString(self: *Tokenizer, quote: u8) Token {
        self.pos += 1; // skip opening quote
        self.col += 1;

        const start = self.pos;

        while (self.pos < self.source.len) {
            const c = self.source[self.pos];

            if (c == quote) {
                // Check if quote is at end of token (followed by whitespace or EOF)
                const next_pos = self.pos + 1;
                if (next_pos >= self.source.len or isWhitespace(self.source[next_pos])) {
                    const value = self.source[start..self.pos];
                    self.pos += 1; // skip closing quote
                    self.col += 1;
                    return .{ .value = value };
                }
            }

            if (c == '\n') {
                self.line += 1;
                self.col = 1;
            } else {
                self.col += 1;
            }
            self.pos += 1;
        }

        // Unterminated string
        return .{ .err = "Unterminated quoted string" };
    }

    /// Read semicolon-delimited text block
    fn readSemicolonText(self: *Tokenizer) Token {
        self.pos += 1; // skip opening semicolon
        self.col += 1;

        // Skip to end of first line
        while (self.pos < self.source.len and self.source[self.pos] != '\n') {
            self.pos += 1;
        }
        if (self.pos < self.source.len) {
            self.pos += 1; // skip newline
            self.line += 1;
            self.col = 1;
        }

        const start = self.pos;

        // Find closing semicolon at beginning of line
        while (self.pos < self.source.len) {
            if (self.col == 1 and self.source[self.pos] == ';') {
                const value = self.source[start..self.pos];
                self.pos += 1; // skip closing semicolon
                self.col += 1;
                return .{ .value = value };
            }

            if (self.source[self.pos] == '\n') {
                self.line += 1;
                self.col = 1;
            } else {
                self.col += 1;
            }
            self.pos += 1;
        }

        // Unterminated text block
        return .{ .err = "Unterminated semicolon text block" };
    }

    /// Read an unquoted value
    fn readUnquotedValue(self: *Tokenizer) Token {
        const start = self.pos;

        while (self.pos < self.source.len) {
            const c = self.source[self.pos];
            if (isWhitespace(c)) break;
            self.pos += 1;
            self.col += 1;
        }

        return .{ .value = self.source[start..self.pos] };
    }

    /// Check if character is whitespace
    fn isWhitespace(c: u8) bool {
        return c == ' ' or c == '\t' or c == '\r' or c == '\n';
    }

    /// Get current line number (1-based)
    pub fn getLine(self: *const Tokenizer) usize {
        return self.line;
    }

    /// Get current column number (1-based)
    pub fn getColumn(self: *const Tokenizer) usize {
        return self.col;
    }
};

/// Check if a value is the CIF "inapplicable" marker
pub fn isInapplicable(value: []const u8) bool {
    return std.mem.eql(u8, value, ".");
}

/// Check if a value is the CIF "unknown" marker
pub fn isUnknown(value: []const u8) bool {
    return std.mem.eql(u8, value, "?");
}

/// Check if a value is null (either . or ?)
pub fn isNull(value: []const u8) bool {
    return isInapplicable(value) or isUnknown(value);
}

// ============================================================================
// Tests
// ============================================================================

test "tokenize data block" {
    var tok = Tokenizer.init("data_1ABC");
    const token = tok.next();
    try std.testing.expectEqual(Token.data_block, std.meta.activeTag(token));
    try std.testing.expectEqualStrings("1ABC", token.data_block);
    try std.testing.expectEqual(Token.eof, std.meta.activeTag(tok.next()));
}

test "tokenize loop" {
    var tok = Tokenizer.init("loop_");
    try std.testing.expectEqual(Token.loop, tok.next());
    try std.testing.expectEqual(Token.eof, std.meta.activeTag(tok.next()));
}

test "tokenize tag" {
    var tok = Tokenizer.init("_atom_site.Cartn_x");
    const token = tok.next();
    try std.testing.expectEqual(Token.tag, std.meta.activeTag(token));
    try std.testing.expectEqualStrings("_atom_site.Cartn_x", token.tag);
}

test "tokenize unquoted value" {
    var tok = Tokenizer.init("123.456");
    const token = tok.next();
    try std.testing.expectEqual(Token.value, std.meta.activeTag(token));
    try std.testing.expectEqualStrings("123.456", token.value);
}

test "tokenize single quoted string" {
    var tok = Tokenizer.init("'hello world'");
    const token = tok.next();
    try std.testing.expectEqual(Token.value, std.meta.activeTag(token));
    try std.testing.expectEqualStrings("hello world", token.value);
}

test "tokenize double quoted string" {
    var tok = Tokenizer.init("\"test value\"");
    const token = tok.next();
    try std.testing.expectEqual(Token.value, std.meta.activeTag(token));
    try std.testing.expectEqualStrings("test value", token.value);
}

test "tokenize semicolon text" {
    const source =
        \\;first line
        \\second line
        \\;
    ;
    var tok = Tokenizer.init(source);
    const token = tok.next();
    try std.testing.expectEqual(Token.value, std.meta.activeTag(token));
    try std.testing.expectEqualStrings("second line\n", token.value);
}

test "skip comments" {
    var tok = Tokenizer.init("# comment\ndata_TEST");
    const token = tok.next();
    try std.testing.expectEqual(Token.data_block, std.meta.activeTag(token));
    try std.testing.expectEqualStrings("TEST", token.data_block);
}

test "tokenize multiple tokens" {
    var tok = Tokenizer.init(
        \\data_TEST
        \\loop_
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\1 C
        \\2 N
    );

    // data_TEST
    var token = tok.next();
    try std.testing.expectEqual(Token.data_block, std.meta.activeTag(token));
    try std.testing.expectEqualStrings("TEST", token.data_block);

    // loop_
    try std.testing.expectEqual(Token.loop, tok.next());

    // _atom_site.id
    token = tok.next();
    try std.testing.expectEqual(Token.tag, std.meta.activeTag(token));
    try std.testing.expectEqualStrings("_atom_site.id", token.tag);

    // _atom_site.type_symbol
    token = tok.next();
    try std.testing.expectEqual(Token.tag, std.meta.activeTag(token));
    try std.testing.expectEqualStrings("_atom_site.type_symbol", token.tag);

    // Values: 1 C 2 N
    token = tok.next();
    try std.testing.expectEqualStrings("1", token.value);
    token = tok.next();
    try std.testing.expectEqualStrings("C", token.value);
    token = tok.next();
    try std.testing.expectEqualStrings("2", token.value);
    token = tok.next();
    try std.testing.expectEqualStrings("N", token.value);

    try std.testing.expectEqual(Token.eof, std.meta.activeTag(tok.next()));
}

test "special values" {
    try std.testing.expect(isInapplicable("."));
    try std.testing.expect(!isInapplicable("?"));
    try std.testing.expect(isUnknown("?"));
    try std.testing.expect(!isUnknown("."));
    try std.testing.expect(isNull("."));
    try std.testing.expect(isNull("?"));
    try std.testing.expect(!isNull("value"));
}

test "case insensitive keywords" {
    var tok = Tokenizer.init("DATA_test LOOP_");
    const token = tok.next();
    try std.testing.expectEqual(Token.data_block, std.meta.activeTag(token));
    try std.testing.expectEqualStrings("test", token.data_block);
    try std.testing.expectEqual(Token.loop, tok.next());
}

test "quoted string with embedded quote char" {
    // In CIF, quotes inside quoted strings work if not followed by whitespace
    var tok = Tokenizer.init("'it''s ok'");
    const token = tok.next();
    try std.testing.expectEqual(Token.value, std.meta.activeTag(token));
    // The embedded '' is kept as-is
    try std.testing.expectEqualStrings("it''s ok", token.value);
}

test "fuzz tokenizer" {
    try std.testing.fuzz({}, struct {
        fn testOne(_: void, input: []const u8) !void {
            var tok = Tokenizer.init(input);
            // Consume all tokens until EOF
            while (true) {
                const token = tok.next();
                if (token == .eof) break;
            }
        }
    }.testOne, .{
        .corpus = &.{
            "data_1ABC",
            "loop_\n_atom_site.id\n_atom_site.type_symbol\n1 C\n2 N",
            "'hello world'",
            "\"test value\"",
            ";first line\nsecond line\n;",
            "# comment\ndata_TEST",
            "DATA_test LOOP_",
            "'it''s ok'",
        },
    });
}
