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

    /// commit will compute the merkle root of the given array.
    /// This root is the vector-commitment.
    pub fn commit (self: *Self, data_ary: [] const [] const u8) !?[hasher.digest_length]u8 {
        if (data_ary.len == 0) {
            return null;
        }

        self.values = try self.allocator.alloc([hasher.digest_length]u8, data_ary.len);

        // Hash the provided values and store them in self.values
        for (data_ary) |datum, idx| {
            var hashed_datum: [hasher.digest_length]u8 = undefined;
            hasher.hash(datum, &hashed_datum, .{});
            self.values.?[idx] = hashed_datum; 
        }
        
        // can be divided by 2
        var pad = try self.allocator.alloc([hasher.digest_length]u8, self.values.?.len);
        defer self.allocator.free(pad);

        // copy over the values. use std.mem.copy instead?
        std.mem.copy([hasher.digest_length]u8, pad, self.values.?);

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

    /// open will return the authenticated-path of a given leaf of the merkle tree
    /// at target_idx
    ///
    /// the caller is responsible for freeing the path returned
    pub fn open(self: *Self, target_idx: usize) ![][hasher.digest_length]u8 {


        if (target_idx >= self.values.?.len) {
            // TODO return error
        }

        var pad = try self.allocator.alloc([hasher.digest_length]u8, self.values.?.len);
        defer self.allocator.free(pad);

        std.mem.copy([hasher.digest_length]u8, pad, self.values.?);

        var path = try self.allocator.alloc(
            [hasher.digest_length]u8, 
            std.math.log2_int_ceil(usize, self.values.?.len) + 1,
        );
        var path_idx: usize = 0;

        // loop over the pad. Each loop represents processing a level
        // in the binary tree. Stop when the last-level of len == 1 is written to.
        var cur_idx = target_idx;
        var len = pad.len;        // length of the level

        std.debug.print("len {}, cur_idx {}\n", .{len, cur_idx});
        while (len > 1) {
            //
            // First, write the authenticating-neighbour into the path
            //

            // to_path_idx is the index of the element in pad 
            // which is the authenticating-neighbour of the target-branch's node.
            var to_path_idx = if (cur_idx % 2 == 0) cur_idx + 1 else cur_idx - 1;

            // if to_path_idx is beyond the length-of-the-level, then we may need 
            // to clone the last element (else-case). This happens when it is not
            // a balanced tree.
            path[path_idx] = if (to_path_idx < len)
                    pad[to_path_idx]
                else
                    pad[to_path_idx - 1];
            path_idx += 1;
            // divide cur_idx by 2 to map it to the next level
            cur_idx = cur_idx >> 1; 

            //
            // Second, write to pad the next level of the binary tree
            //
            var read_idx: usize = 0; 
            var write_idx: usize = 0; 
            while (read_idx < len) {
                pad[write_idx] = hashAdjacent(pad, read_idx, len);
                read_idx += 2;
                write_idx += 1;
            }
            len = write_idx;
        }

        
        path = self.allocator.shrink(path, path_idx);
        return path;
    }

    /// verify accepts the target value, target_idx and the path. Using these, it will
    /// walk up the branch from the target to the root and re-calculate the root.
    ///
    /// Finally, it will verify that its calculated root is the same as the root passed in.
    ///
    /// Why? This verifies that a given leaf (target value) is an element of the
    /// committed vector (root) at the given index (target_idx)
    pub fn verify(
        self: *Self, 
        root: [hasher.digest_length]u8,
        _target_idx: usize,
        path: [][hasher.digest_length]u8,
        target: []const u8,
    ) !bool {
        var target_idx = _target_idx;

        if (self.values.?.len <= 0) {
            return false;
        }
        
        var pad = try self.allocator.alloc([hasher.digest_length]u8, self.values.?.len);
        defer self.allocator.free(pad);

        std.mem.copy([hasher.digest_length]u8, pad, self.values.?);
        
        var target_hash: [hasher.digest_length]u8 = undefined;
        hasher.hash(target, &target_hash, .{});

        var path_idx: usize = 0;
        while (path_idx < path.len) {
            var path_val = path[path_idx];
    
            var buf = if (target_idx % 2 == 0) 
                    copyToBuf(target_hash, path_val) 
                else 
                    copyToBuf(path_val, target_hash);
            hasher.hash(buf[0..], &target_hash, .{});

            target_idx = target_idx >> 1; // divide by 2
            path_idx += 1;
        }

        return std.mem.eql(u8, root[0..], target_hash[0..]);
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

fn doTest(data_ary: []const []const u8, idx: usize, path_len: usize) anyerror!void { 
    var m = SlimMerkleTrie.init(testing.allocator);
    defer m.deinit();

    const root = try m.commit(data_ary);
    try std.testing.expect(root != null);

    var path = try m.open(idx);
    defer testing.allocator.free(path);
    try std.testing.expect(path.len == path_len);

    var is_valid = try m.verify(root.?, idx, path, data_ary[idx]);
    try std.testing.expect(is_valid);
}