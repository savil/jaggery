const std = @import("std");
const testing = std.testing;
const Allocator = std.mem.Allocator;
// TODO savil. Should this import the library instead of the file?
// TODO savil. Rename FieldElement.zig to FiniteFields.zig?
const FiniteFields = @import("finite-fields");
const FieldElement = FiniteFields.FieldElement;
const Field = FiniteFields.Field;


const Polynomial = struct {
    const Self = @This();

    allocator: *const Allocator,
    coefficients: []FieldElement,
    
    // TODO savil: we need an allocator to store the coefficients in the heap
    pub fn init(allocator: *const Allocator, coefficients: []const FieldElement) !Self {
        // make internal copy of coefficients since we own it now
        var coeffs = try allocator.alloc(FieldElement, coefficients.len);
        std.mem.copy(FieldElement, coeffs[0..], coefficients[0..]);

        return Self {
            .allocator = allocator,
            .coefficients = coeffs,
        };
    }

    pub fn deinit(self: Self) void {
        self.allocator.free(self.coefficients);
    }

    fn degree(self: Self) ?usize {
        if (self.coefficients.len == 0) {
            return null;
        }

        var all_zeros = true;
        for (self.coefficients) | coeff | {
            if (!coeff.isZero()) {
                all_zeros = false;
            }
        }
        if (all_zeros) {
            return null;
        }

        var maxIndex: usize = 0;
        for (self.coefficients) | coeff, idx | {
            if (!coeff.isZero()) {
                maxIndex = idx;
            }
        }
        return maxIndex; 
    }

    fn neg(self: *Self) *Self {
        for (self.coefficients) | coeff, idx | {
            self.coefficients[idx] = coeff.neg();
        }
        return self;
    }

    // add appends the coefficients of other onto the coefficients of add.
    fn add(self: *Self, other: *const Polynomial) !*Self {

        if (other.degree() == null) {
            return self;
        }

        const combined_len = std.math.max(self.coefficients.len, other.coefficients.len);
        const coeffs = try self.allocator.alloc(FieldElement, combined_len);
        std.mem.copy(FieldElement, coeffs[0..], self.coefficients);
        for (other.coefficients) | coeff, idx | {
            if (idx < self.coefficients.len) {
                coeffs[idx] = coeffs[idx].add(coeff);
            } else {
                coeffs[idx] = coeff;
            }
        }
        self.allocator.free(self.coefficients);
        self.coefficients = coeffs;
        return self;
    }

    fn sub(self: *Self, other: *const Polynomial) !*Self {
        var other2 = try Polynomial.init(self.allocator, other.coefficients);
        defer other2.deinit();
        return try self.add(other2.neg());
    }

    fn mul(self: *Self, other: *const Polynomial) !*Self {
        if ((self.coefficients.len == 0) or (other.coefficients.len == 0)) {
            var buf = try self.allocator.alloc(FieldElement, 0);
            self.allocator.free(self.coefficients);
            self.coefficients = buf;
            return self;
        }

        const buf_size = self.coefficients.len + other.coefficients.len - 1;
        var buf = try self.allocator.alloc(FieldElement, buf_size);
        for (buf) | _, idx | {
            buf[idx] = self.coefficients[0].field.zero();
        }

        for (self.coefficients) | self_coeff, i | {
            for (other.coefficients) | other_coeff, j | {
                buf[i + j] = buf[i + j].add(self_coeff.mul(other_coeff));
            }
        }

        self.allocator.free(self.coefficients);
        self.coefficients = buf;
        return self;
    }

    fn eq(self: *Self, other: *const Polynomial) bool {
        const self_degree = self.degree();
        const other_degree = other.degree();
        if (self_degree == null and other_degree == null) {
            return true;
        }
        if (self_degree == null or other_degree == null) {
            return false;
        }
        if (self_degree.? != other_degree.?) {
            return false;
        }

        for (self.coefficients) | self_coeff, idx | {
            if (!self_coeff.eq(other.coefficients[idx])) {
                return false;
            }
        }
        return true;
    }

    fn neq(self: *Self, other: *const Polynomial) bool {
        return !self.eq(other);
    }

    fn isZero(self: *Self) bool {
        return self.degree() == null;
    }

    fn leadingCoefficient(self: *Self) ?FieldElement {
        const d = self.degree();
        if (d == null) {
            return null;
        }
        return self.coefficients[d.?];
    }
};

test "init polynomial" {
    const field = Field.init(19);
    const fe = FieldElement.init(2, field);
    const fes = [_]FieldElement{fe};
    const polynomial = try Polynomial.init(&testing.allocator, fes[0..]);
    defer polynomial.deinit();

    _ = polynomial;
    try testing.expect(true); // just want the above to compile.
}

test "degrees-zero-len" {
    const fes = [_]FieldElement{};
    const polynomial = try Polynomial.init(&testing.allocator, fes[0..]);
    defer polynomial.deinit();
    try testing.expect(polynomial.degree() == null);
}

test "degrees-all-zeros" {
    const field = Field.init(19);
    const fe1 = FieldElement.init(0, field);
    const fe2 = FieldElement.init(0, field);
    const fes = [_]FieldElement{fe1, fe2};
    const polynomial = try Polynomial.init(&testing.allocator, fes[0..]);
    defer polynomial.deinit();
    try testing.expect(polynomial.degree() == null);
}

test "degrees-max-index" {
    const field = Field.init(19);
    const fe1 = FieldElement.init(0, field);
    const fe2 = FieldElement.init(3, field);
    const fe3 = FieldElement.init(0, field);
    const fe4 = FieldElement.init(5, field);
    const fe5 = FieldElement.init(0, field);
    const fes = [_]FieldElement{fe1, fe2, fe3, fe4, fe5};
    const polynomial = try Polynomial.init(&testing.allocator, fes[0..]);
    defer polynomial.deinit();
    try testing.expect(polynomial.degree().? == 3);
}

test "neg" {
    const field = Field.init(19);
    const fe1 = FieldElement.init(0, field);
    const fe2 = FieldElement.init(3, field);
    const fes = [_]FieldElement{fe1, fe2};
    var polynomial = try Polynomial.init(&testing.allocator, fes[0..]);
    defer polynomial.deinit();

    _ = polynomial.neg();
    try testing.expect(polynomial.degree().? == 1);
    try testing.expect(polynomial.coefficients[1].eq(fe2.neg()));
}

test "add empty" {
    const field = Field.init(19);
    const fe1 = FieldElement.init(0, field);
    const fe2 = FieldElement.init(3, field);
    const fes = [_]FieldElement{fe1, fe2};
    var polynomial1 = try Polynomial.init(&testing.allocator, fes[0..]);
    defer polynomial1.deinit();

    try testing.expect(polynomial1.degree().? == 1);

    var polynomial2 = try Polynomial.init(&testing.allocator, &[_]FieldElement{});
    defer polynomial2.deinit();
    try testing.expect(polynomial2.degree() == null);

    // add empty-polynomial to polynomial1. Ensure polynomial1 is same.
    _ = try polynomial1.add(&polynomial2);
    try testing.expect(polynomial1.degree().? == 1);
    try testing.expect(polynomial1.coefficients[0].eq(fe1));
    try testing.expect(polynomial1.coefficients[1].eq(fe2));

    // add polynomial1 to the empty-polynomial. Now, polynomial2 is no longer empty.
    _ = try polynomial2.add(&polynomial1);
    try testing.expect(polynomial2.degree().? == 1);
    try testing.expect(polynomial2.coefficients[0].eq(fe1));
    try testing.expect(polynomial2.coefficients[1].eq(fe2));
}

test "add" {
    const field = Field.init(19);
    const fe1 = FieldElement.init(0, field);
    const fe2 = FieldElement.init(3, field);
    const fes = [_]FieldElement{fe1, fe2, fe2};
    var polynomial1 = try Polynomial.init(&testing.allocator, fes[0..]);
    defer polynomial1.deinit();

    try testing.expect(polynomial1.degree().? == 2);

    const polynomial2 = try Polynomial.init(&testing.allocator, &[_]FieldElement{fe2, fe1, fe2});
    defer polynomial2.deinit();
    try testing.expect(polynomial2.degree().? == 2);
    
    _ = try polynomial1.add(&polynomial2);
    try testing.expect(polynomial1.degree().? == 2);
    try testing.expect(polynomial1.coefficients[0].eq(fe2));
    try testing.expect(polynomial1.coefficients[1].eq(fe2));
    try testing.expect(polynomial1.coefficients[2].eq(FieldElement.init(6, field)));
}

test "sub" {
    const field = Field.init(19);
    const fe1 = FieldElement.init(0, field);
    const fe2 = FieldElement.init(3, field);
    const fes = [_]FieldElement{fe1, fe2, fe2};
    var polynomial1 = try Polynomial.init(&testing.allocator, fes[0..]);
    defer polynomial1.deinit();

    try testing.expect(polynomial1.degree().? == 2);

    const polynomial2 = try Polynomial.init(&testing.allocator, &[_]FieldElement{fe2, fe1, fe2});
    defer polynomial2.deinit();
    try testing.expect(polynomial2.degree().? == 2);
    
    _ = try polynomial1.sub(&polynomial2);
    try testing.expect(polynomial1.degree().? == 1);
}

test "mul" {

    const field = Field.init(19);
    const fe1 = FieldElement.init(0, field);
    const fe2 = FieldElement.init(3, field);
    const fes = [_]FieldElement{fe1, fe2, fe2};
    var polynomial1 = try Polynomial.init(&testing.allocator, fes[0..]);
    defer polynomial1.deinit();

    try testing.expect(polynomial1.degree().? == 2);

    const polynomial2 = try Polynomial.init(&testing.allocator, &[_]FieldElement{fe2, fe1, fe2});
    defer polynomial2.deinit();
    try testing.expect(polynomial2.degree().? == 2);

    _ = try polynomial1.mul(&polynomial2);
    try testing.expect(polynomial1.degree().? == 4);
}

test "mul empty" {

    const field = Field.init(19);
    const fe1 = FieldElement.init(0, field);
    const fe2 = FieldElement.init(3, field);
    const fes = [_]FieldElement{fe1, fe2, fe2};
    var polynomial1 = try Polynomial.init(&testing.allocator, fes[0..]);
    defer polynomial1.deinit();

    try testing.expect(polynomial1.degree().? == 2);

    var polynomial2 = try Polynomial.init(&testing.allocator, &[_]FieldElement{});
    defer polynomial2.deinit();
    try testing.expect(polynomial2.degree() == null);
    
    _ = try polynomial2.mul(&polynomial1);
    try testing.expect(polynomial2.degree() == null);
    try testing.expect(polynomial1.degree().? == 2);

    _ = try polynomial1.mul(&polynomial2);
    try testing.expect(polynomial1.degree() == null);
    try testing.expect(polynomial2.degree() == null);
}

test "eq - both empty" {
    // First, make both polynomials empty
    var polynomial1 = try Polynomial.init(&testing.allocator, &[_]FieldElement{});
    defer polynomial1.deinit();
    try testing.expect(polynomial1.degree() == null);

    var polynomial2 = try Polynomial.init(&testing.allocator, &[_]FieldElement{});
    defer polynomial2.deinit();
    try testing.expect(polynomial2.degree() == null);

    try testing.expect(polynomial1.eq(&polynomial2));
}

test "eq - one empty" {
    const field = Field.init(19);
    const fe1 = FieldElement.init(0, field);
    const fe2 = FieldElement.init(3, field);
    const fes = [_]FieldElement{fe1, fe2, fe2};
    var polynomial1 = try Polynomial.init(&testing.allocator, fes[0..]);
    defer polynomial1.deinit();

    try testing.expect(polynomial1.degree().? == 2);

    var polynomial2 = try Polynomial.init(&testing.allocator, &[_]FieldElement{});
    defer polynomial2.deinit();
    try testing.expect(polynomial2.degree() == null);

    try testing.expect(polynomial1.neq(&polynomial2));

    try testing.expect(polynomial2.neq(&polynomial1));
}

test "eq - both equal" {

    const field = Field.init(19);
    const fe1 = FieldElement.init(0, field);
    const fe2 = FieldElement.init(3, field);
    const fes = [_]FieldElement{fe1, fe2, fe2};
    var polynomial1 = try Polynomial.init(&testing.allocator, fes[0..]);
    defer polynomial1.deinit();

    var polynomial2 = try Polynomial.init(&testing.allocator, fes[0..]);
    defer polynomial2.deinit();
    try testing.expect(polynomial1.eq(&polynomial2));
}

test "isZero" {
    const field = Field.init(19);
    const fe1 = FieldElement.init(0, field);
    const fe2 = FieldElement.init(3, field);
    const fes = [_]FieldElement{fe1, fe2, fe2};
    var polynomial1 = try Polynomial.init(&testing.allocator, fes[0..]);
    defer polynomial1.deinit();
    try testing.expect(!polynomial1.isZero());

    var polynomial2 = try Polynomial.init(&testing.allocator, &[_]FieldElement{});
    defer polynomial2.deinit();
    try testing.expect(polynomial2.isZero());
}

test "leadingCoefficient" {
    const field = Field.init(19);
    const fe1 = FieldElement.init(0, field);
    const fe2 = FieldElement.init(3, field);
    const fes = [_]FieldElement{fe1, fe2, fe1};
    var polynomial1 = try Polynomial.init(&testing.allocator, fes[0..]);
    defer polynomial1.deinit();
    try testing.expect(!polynomial1.isZero());
    try testing.expect(polynomial1.leadingCoefficient().?.eq(fe2));

    var polynomial2 = try Polynomial.init(&testing.allocator, &[_]FieldElement{});
    defer polynomial2.deinit();
    try testing.expect(polynomial2.isZero());
    try testing.expect(polynomial2.leadingCoefficient() == null);
}