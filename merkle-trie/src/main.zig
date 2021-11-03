const std = @import("std");
const testing = std.testing;
const Allocator = std.mem.Allocator;

/// fun inspiration:
/// https://aszepieniec.github.io/stark-anatomy/basic-tools#merkle-tree
const merkleTrie = struct {

    const Self = @This();
    const hasher = std.crypto.hash.blake2.Blake2b256;

    values: std.ArrayList([hasher.digest_length]u8),
    num_items: u64,

    pub fn init(allocator: *Allocator) Self {
        return Self{
            .values = std.ArrayList([hasher.digest_length]u8).init(allocator),
            .num_items = 0,
        };
    }

    pub fn deinit(self: Self) void {
        self.values.deinit();
    }

    /// commit eagerly builds a balanced binary trie. This is not space-efficient
    /// but is simpler for me to reason about :-). Optimizations can be made later
    /// once this works.
    pub fn commit(self: *Self, dataAry: [] const []const u8) !void {
        self.num_items = dataAry.len;

        for (dataAry) |datum| {
            var hashed_datum: [hasher.digest_length]u8 = undefined;
            hasher.hash(datum, &hashed_datum, .{});
            try self.values.append(hashed_datum); 
        }
        try self.ensureValuesLengthIsPowerOfTwo();


        var level_start_idx: usize = 0;
        var level_idx: usize = 0;
        var level_len: usize = self.values.items.len;
        while (level_len >= 2) {
            while (level_idx < level_len) {
                var current_idx = level_start_idx + level_idx;
                var buf = self.copyToBuf(current_idx);
                var hashed_datum: [hasher.digest_length]u8 = undefined;
                hasher.hash(buf[0..], &hashed_datum, .{});
                try self.values.append(hashed_datum); 

                // copyToBuf processes two array elements, 
                // so we move ahead by two.
                level_idx += 2; 
            }

            level_start_idx += level_len;
            level_idx = 0;
            level_len = level_len / 2;
        }
    }

    pub fn open(idx: u64) void {
        _ = idx;

    }

    pub fn verify() void {

    }

    fn copyToBuf(self: *Self, cp_idx: usize) []const u8 {

        var buf: [2* hasher.digest_length] u8 = undefined;
        var idx: usize = 0;
        for (self.values.items[cp_idx]) | char | {
            buf[idx] = char;
            idx += 1;
        }
        for (self.values.items[cp_idx + 1]) | char | {
            buf[idx] = char;
            idx += 1;
        }
        return buf[0..];
    }

    fn ensureValuesLengthIsPowerOfTwo(self: *Self) !void {
        while (!self.valuesLengthIsPowerOfTwo()) {
            const datum = "";
            var hashed_datum: [hasher.digest_length]u8 = undefined;
            hasher.hash(datum, &hashed_datum, .{});
            try self.values.append(hashed_datum);
        }
    }

    // https://bits.stephan-brumme.com/isPowerOfTwo.html
    fn valuesLengthIsPowerOfTwo(self: *Self) bool {
        const lenVals = self.values.items.len;
        if (lenVals < 2) {
            return true;
        }
        return (lenVals & (lenVals - 1)) == 0;
    }
};

test "init empty merkle trie" {
    const m = merkleTrie.init(testing.allocator);
    defer m.deinit();
}

test "empty merkle trie" {
    var m = merkleTrie.init(testing.allocator);
    defer m.deinit();
    
    const dataAry = [_][]const u8{};
    try m.commit(dataAry[0..]);
}

test "(small) merkle trie" {
    var m = merkleTrie.init(testing.allocator);
    defer m.deinit();

    const dataAry = [4][]const u8{"savil", "diana", "noorjahan", "ruksaar"};
    try m.commit(dataAry[0..]);
}