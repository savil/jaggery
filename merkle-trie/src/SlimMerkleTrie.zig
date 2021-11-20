const std = @import("std");
const testing = std.testing;
const Allocator = std.mem.Allocator;


const SlimMerkleTrie = struct {
    const Self = @This();

    const hasher = std.crypto.hash.blake2.Blake2b256;

    allocator: *Allocator,
    values: [][hasher.digest_length]u8,

    pub fn init(allocator: *Allocator) Self {
        return Self {
            .allocator = allocator,
            .values = undefined,
        };
    }

    pub fn deinit(self: Self) void {
        _ = self;
        // do equivalent of self.allocator.deinit(); ?
    }   

    pub fn commit (self: *Self, data_ary: [] const [] const u8) !?[hasher.digest_length]u8 {
        self.values = try self.allocator.alloc(
            [hasher.digest_length]u8, 
            hasher.digest_length * data_ary.len,
        );

        return null;
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