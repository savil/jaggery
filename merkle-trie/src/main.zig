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

    values: std.ArrayList([]u8),

    pub fn init(allocator: *Allocator) Self {
        const list = std.ArrayList([]u8).init(allocator);
        return Self{
            .values = list
        };
    }

    pub fn deinit(self: Self) void {
        self.values.deinit();
    }

    pub fn commit() void {

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