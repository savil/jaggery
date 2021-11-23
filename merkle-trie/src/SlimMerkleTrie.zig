const std = @import("std");
const testing = std.testing;
const Allocator = std.mem.Allocator;


const SlimMerkleTrie = struct {
    const Self = @This();

    const hasher = std.crypto.hash.blake2.Blake2b256;

    allocator: *Allocator,
    values: ?[][hasher.digest_length]u8,

    pub fn init(allocator: *Allocator) Self {
        return Self {
            .allocator = allocator,
            .values = null,
        };
    }

    pub fn deinit(self: Self) void {
        if (self.values != null) {
            self.allocator.free(self.values.?);
        }
    }   

    pub fn commit (self: *Self, data_ary: [] const [] const u8) !?[hasher.digest_length]u8 {
        self.values = try self.allocator.alloc(
            [hasher.digest_length]u8, 
            hasher.digest_length * data_ary.len,
        );

        // Hash the provided values and store them in self.values
        for (data_ary) |datum, idx| {
            var hashed_datum: [hasher.digest_length]u8 = undefined;
            hasher.hash(datum, &hashed_datum, .{});
            self.values.?[idx] = hashed_datum; 
        }

        var pad = try self.allocator.alloc(
            [hasher.digest_length]u8,
            hasher.digest_length * data_ary.len, // can be divided by 2
        );
        defer self.allocator.free(pad);

        // copy over the values. use std.mem.copy instead?
        for (self.values.?) | value, idx | {
            pad[idx] = value;
        }

        if (pad.len <= 0) {
            return null;
        }

        // loop over the pad. Each loop represents processing a level
        // in the binary tree. Stop when the last-level of len == 1 is written to.
        var len = pad.len;        // length of the level
        while (len > 1) {
            var read_idx: usize = 0; 
            var write_idx: usize = 0; 

            while (read_idx < len) {
                pad[write_idx] = hashAdjacent(pad, read_idx, len);
                read_idx += 2;
                write_idx += 1;
            }

            len = write_idx;
        }
        return pad[0];
    }

    pub fn open(self: *Self, idx: usize) ![][hasher.digest_length]u8 {
        _ = self;
        _ = idx;
        
        return undefined;
    }

    pub fn verify(
        self: *Self, 
        root: [hasher.digest_length]u8,
        idx: usize,
        path: [][hasher.digest_length]u8,
        target: []const u8,
    ) bool {
        _ = self;
        _ = root;
        _ = idx;
        _ = path;
        _ = target;
        return false;
    }

    fn hashAdjacent(pad: [][hasher.digest_length]u8, read_idx: usize, len: usize) [hasher.digest_length]u8 {
        // if one element is dangling, then we hash it to itself.
        var read_adj_idx = read_idx + 1;
        if (read_adj_idx >= len) {
            read_adj_idx = read_idx;
        }

        // Hash each pair of adjacent elements
        var buf = copyToBuf(pad[read_idx], pad[read_adj_idx]);
        var hashed_datum: [hasher.digest_length]u8 = undefined;
        hasher.hash(buf[0..], &hashed_datum, .{});
        return hashed_datum;
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
};


test "empty merkle trie" {
    var m = SlimMerkleTrie.init(testing.allocator);
    defer m.deinit();
    
    const data_ary = [_][]const u8{};
    _ = try m.commit(data_ary[0..]);
}

test "size-one merkle trie" {
    const data_ary = [_][]const u8{"asterix"};
    try doTest(data_ary[0..], 0, 0);
}

fn doTest(data_ary: []const []const u8, idx: usize, path_len: usize) anyerror!void { 
    var m = SlimMerkleTrie.init(testing.allocator);
    defer m.deinit();

    const root = try m.commit(data_ary);
    try std.testing.expect(root != null);

    var path = try m.open(idx);
    try std.testing.expect(path.len == path_len);

    var is_valid = m.verify(root.?, idx, path, data_ary[idx]);
    try std.testing.expect(is_valid);
}