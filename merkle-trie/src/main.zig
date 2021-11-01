const std = @import("std");
const testing = std.testing;
const Allocator = std.mem.Allocator;

//fn newMerkleTrie() merkleTrie {
 //   return merkleTrie{
 //       .values = nil
  //  };
//}

/// fun inspiration:
/// https://aszepieniec.github.io/stark-anatomy/basic-tools#merkle-tree
const merkleTrie = struct {

    const Self = @This();
    const hasher = std.crypto.hash.blake2.Blake2b256;

    values: std.ArrayList([hasher.digest_length]u8),

    pub fn init(allocator: *Allocator) Self {
        return Self{
            .values = std.ArrayList([hasher.digest_length]u8).init(allocator),
        };
    }

    pub fn deinit(self: Self) void {
        self.values.deinit();
    }

    pub fn commit(self: *Self, dataAry: [] const []const u8) !void {
        _ = self;
        _ = dataAry;

        for (dataAry) |datum| {
            var hashed_datum: [hasher.digest_length]u8 = undefined;
            hasher.hash(datum, &hashed_datum, .{});
            try self.values.append(hashed_datum); 
        }
    }

    pub fn open() void {

    }

    pub fn verify() void {

    }

};

test "init empty merkle trie" {
    const m = merkleTrie.init(testing.allocator);
    defer m.deinit();
}

test "commit merkle trie" {
    var m = merkleTrie.init(testing.allocator);
    defer m.deinit();

    const dataAry = [2][]const u8{"savil", "diana"};
    try m.commit(dataAry[0..]);
}