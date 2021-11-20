
const std = @import("std");
const testing = std.testing;
const Allocator = std.mem.Allocator;


/// fun inspiration:
/// https://aszepieniec.github.io/stark-anatomy/basic-tools#merkle-tree
///
/// This implementation is not memory efficient. In fact, it grows exponentially
/// because we force the size of the internal ArrayList to be a power of two! Yikes.
///
/// Maybe this would be a good way to explore how "interfaces" may work in ziglang. 
/// I could keep this inefficient implementation and have a new struct SmallerMemoryMerkleTrie
/// that I also test against.
const SpaceyMerkleTrie = struct {

    const Self = @This();
    const hasher = std.crypto.hash.blake2.Blake2b256;

    allocator: *Allocator,
    values: std.ArrayList([hasher.digest_length]u8),
    num_items: u64, // total number of items added initially to the tree

    pub fn init(allocator: *Allocator) Self {
        return Self{
            .allocator = allocator,
            .values = std.ArrayList([hasher.digest_length]u8).init(allocator),
            .num_items = 0,
        };
    }

    pub fn deinit(self: Self) void {
        self.values.deinit();
    }

    /// commit eagerly builds a balanced binary trie. This is not space-efficient
    /// but is simpler for me to reason about :-)
    pub fn commit(self: *Self, data_ary: [] const []const u8) !?[hasher.digest_length]u8 {

        // Hash the provided values and store them in self.values
        for (data_ary) |datum| {
            var hashed_datum: [hasher.digest_length]u8 = undefined;
            hasher.hash(datum, &hashed_datum, .{});
            try self.values.append(hashed_datum); 
        }
        try self.ensureValuesLengthIsPowerOfTwo();
        self.num_items = self.values.items.len;


        //
        // Loop over every level in the binary tree, starting from the bottom level
        // until the root.
        //

        var level_start_idx: usize = 0; 
        var level_idx: usize = 0; // idx within the level, starting at level_start_idx
        var level_len: usize = self.values.items.len;

        // Increment the level. 
        // The last level processed is the one that adds the root-element to self.values
        // Hence we don't process the top-most level (i.e. level_len == 1)
        while (level_len >= 2) {

            // Loop within the level:
            while (level_idx < level_len) {

                // Hash each pair of adjacent elements
                var current_idx = level_start_idx + level_idx;
                var buf = copyToBuf(
                    self.values.items[current_idx], 
                    self.values.items[current_idx + 1],
                );
                var hashed_datum: [hasher.digest_length]u8 = undefined;
                hasher.hash(buf[0..], &hashed_datum, .{});
                try self.values.append(hashed_datum); 

                // copyToBuf processes two array elements, 
                // so we move ahead by two.
                level_idx += 2; 
            }

            // Reset these variables to the next binary-tree level
            level_start_idx += level_len;
            level_idx = 0;
            level_len = level_len / 2;
        }

        // TODO savil. Move this to the top?
        if (self.values.items.len == 0) {
            return null;
        }
        return self.values.items[self.values.items.len - 1];
    }

    // Open up the tree, starting from the `idx` element at the bottom-level.
    pub fn open(self: *Self, idx: u64) !std.ArrayList([hasher.digest_length]u8) {

        // TODO savil. Throw error if idx > length of added-values.
        if (idx >= self.num_items) {
            // TODO return an error
        }

        var path = std.ArrayList([SpaceyMerkleTrie.hasher.digest_length]u8).init(self.allocator);

        var level_start_idx: usize = 0;
        var level_idx: usize = idx;
        var level_len: usize = self.num_items;
        while (level_len > 1) {
            var path_idx = if (level_idx % 2 == 0) level_idx + 1 else level_idx - 1;
            try path.append(self.values.items[path_idx]);
            
            level_start_idx += level_len;
            level_idx = level_start_idx + @divTrunc(level_idx, 2);
            
            level_len = level_len / 2;
        }
        return path;
    }

    pub fn verify(
        self: *Self,
        root: [hasher.digest_length]u8, 
        idx: usize, 
        path: std.ArrayList([hasher.digest_length]u8), 
        target: []const u8,
    ) bool {

        var target_hash: [hasher.digest_length]u8 = undefined;
        hasher.hash(target, &target_hash, .{});

        var path_idx: usize = 0;
        var level_start_idx: usize = 0;
        var level_idx: usize = idx;
        var level_len: usize = self.num_items;
        var buf: []const u8 = "";

        while (level_len > 1) {
            var item = path.items[path_idx];
            buf = if (level_idx % 2 == 0) 
                    copyToBuf(target_hash, item)
                else 
                    copyToBuf(item, target_hash);
                
            hasher.hash(buf[0..], &target_hash, .{});
            

            path_idx += 1;
            level_start_idx += level_len;
            level_idx = level_start_idx + @divTrunc(level_idx, 2);
            level_len = level_len / 2;
        }

        return std.mem.eql(u8, root[0..], target_hash[0..]);
    }

    // lol, surely there is a way to make std.mem.copy do this:
    fn copyToBuf(item1: [hasher.digest_length]u8, item2: [hasher.digest_length]u8) []const u8 {

        var buf: [2* hasher.digest_length] u8 = undefined;
        var idx: usize = 0;
        for (item1) | char | {
            buf[idx] = char;
            idx += 1;
        }
        for (item2) | char | {
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
        const len_vals = self.values.items.len;
        if (len_vals < 2) {
            return true;
        }
        return (len_vals & (len_vals - 1)) == 0;
    }
};

test "init empty merkle trie" {
    const m = SpaceyMerkleTrie.init(testing.allocator);
    defer m.deinit();
}

test "empty merkle trie" {
    var m = SpaceyMerkleTrie.init(testing.allocator);
    defer m.deinit();
    
    const data_ary = [_][]const u8{};
    _ = try m.commit(data_ary[0..]);

    //try std.testing.expect(true == false);
}

test "size-one merkle trie" {
    const data_ary = [_][]const u8{"asterix"};
    try doTest(data_ary[0..], 0, 0);
}

test "size-two merkle trie" {
    const data_ary = [_][]const u8{"asterix", "getafix"};
    try doTest(data_ary[0..], 0, 1);
    try doTest(data_ary[0..], 1, 1);
}

test "size-three merkle trie" {
    const data_ary = [_][]const u8{"asterix", "getafix", "obelix"};
    try doTest(data_ary[0..], 0, 2);
    try doTest(data_ary[0..], 1, 2);
    try doTest(data_ary[0..], 2, 2);
}

test "size-four merkle trie" {
    const data_ary = [_][]const u8{"asterix", "getafix", "obelix", "cacofonix"};
    try doTest(data_ary[0..], 0, 2);
    try doTest(data_ary[0..], 1, 2);
    try doTest(data_ary[0..], 2, 2);
    try doTest(data_ary[0..], 3, 2);
}

fn doTest(data_ary: []const []const u8, idx: u64, path_len: usize) anyerror!void { 
    var m = SpaceyMerkleTrie.init(testing.allocator);
    defer m.deinit();

    const root = try m.commit(data_ary);
    try std.testing.expect(root != null);

    var path = try m.open(idx);
    defer path.deinit();
    try std.testing.expect(path.items.len == path_len);

    var is_valid = m.verify(root.?, idx, path, data_ary[idx]);
    try std.testing.expect(is_valid);
}
