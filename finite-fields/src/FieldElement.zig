const std = @import("std");
const assert = std.debug.assert;

/// This file has FieldElement and Fields.

const FieldElement = struct {
    const Self = @This();

    const TValue = i256;

    value: TValue,
    field: Field,

    pub fn init(value: TValue, field: Field) Self {
        return Self {
            .value = value,
            .field = field,
        };
    }

    pub fn add(self: Self, right: FieldElement) FieldElement {
        return self.field.add(self, right);
    }

    pub fn mul(self: Self, right: FieldElement) FieldElement {
        return self.field.multiply(self, right);
    }    

    pub fn sub(self: Self, right: FieldElement) FieldElement {
        return self.field.subtract(self, right);
    }

    pub fn div(self: Self, right: FieldElement) FieldElement {
        return self.field.divide(self, right);
    }

    pub fn neg(self: Self) FieldElement {
        return self.field.negate(self);
    }

    pub fn inv(self: Self) !FieldElement {
        return self.field.inverse(self);
    }

    // Note, in Alan's implementation this is __xor__ def
    // because C languages use ^ for xor but he probably actually
    // wanted to implement exponentiation (as in the math operator)
    //
    // https://www.khanacademy.org/computing/computer-science/cryptography/modarithmetic/a/fast-modular-exponentiation
    pub fn exp(self: Self, exponent: i64) FieldElement {
        var powOfTwo: i64 = 1;
        var valuePowOfTwo = FieldElement.init(self.value, self.field);
        var result = self.field.one();

        while (powOfTwo <= exponent) {
            if (powOfTwo & exponent == powOfTwo) {
                result = result.mul(valuePowOfTwo);
            }
            valuePowOfTwo = valuePowOfTwo.mul(valuePowOfTwo);

            powOfTwo = powOfTwo << 1;
        }
        return result;
    }

    pub fn eq(self: Self, other: FieldElement) FieldElement {
        return self.value == other.value;
    }

    pub fn neq(self: Self, other: FieldElement) FieldElement {
        return self.value != other.value;
    }

    pub fn bytes() void {
        // TODO
    }

    pub fn isZero(self: Self) bool {
        return self.value == 0;
    }

};


test "init field element" {
    var field = Field.init(3);

    var fieldElem = FieldElement.init(5, field);
    _ = fieldElem;
}

test "add field elements" {
    var field = Field.init(7);

    var fe0 = field.zero();
    var fe1 = field.one();
    var fe3 = FieldElement.init(3, field);
    var fe4 = FieldElement.init(4, field);
    var fe5 = FieldElement.init(5, field);
    var fe7 = FieldElement.init(7, field);

    try std.testing.expectEqual(fe0.add(fe0), fe0);
    try std.testing.expectEqual(fe7.add(fe7), fe0);

    try std.testing.expectEqual(fe1.add(fe3), fe4);

    // (5 + 3) % 7 == 8 % 7 == 1
    try std.testing.expectEqual(fe5.add(fe3), fe1);
}

test "multiply field elements" {
    var field = Field.init(7);

    var fe0 = field.zero();
    var fe1 = field.one();
    var fe3 = FieldElement.init(3, field);
    //var fe4 = FieldElement.init(4, field);
    var fe5 = FieldElement.init(5, field);
    var fe7 = FieldElement.init(7, field);

    try std.testing.expectEqual(fe0.mul(fe0), fe0);
    try std.testing.expectEqual(fe7.mul(fe7), fe0);

    try std.testing.expectEqual(fe1.mul(fe3), fe3);

    // (5 * 3) % 7 == 15 % 7 == 1
    try std.testing.expectEqual(fe5.mul(fe3), fe1);
}

test "subtract field elements" {
    var field = Field.init(7);

    var fe0 = field.zero();
    var fe1 = field.one();
    var fe2 = FieldElement.init(2, field);
    var fe3 = FieldElement.init(3, field);
    var fe4 = FieldElement.init(4, field);
    var fe5 = FieldElement.init(5, field);
    var fe6 = FieldElement.init(6, field);
    var fe7 = FieldElement.init(7, field);

    try std.testing.expectEqual(fe0.sub(fe0), fe0);
    try std.testing.expectEqual(fe2.sub(fe0), fe2);
    try std.testing.expectEqual(fe7.sub(fe0), fe0);

    try std.testing.expectEqual(fe7.sub(fe1), fe6);
    try std.testing.expectEqual(fe7.sub(fe2), fe5);
    try std.testing.expectEqual(fe7.sub(fe3), fe4);

    try std.testing.expectEqual(fe2.sub(fe5), fe4); // (2 - 5) % 7 == -3 % 7 == 4
}

test "negate field elements" {
    var field = Field.init(7);

    var fe0 = field.zero();
    var fe1 = field.one();
    var fe2 = FieldElement.init(2, field);
    var fe3 = FieldElement.init(3, field);
    var fe4 = FieldElement.init(4, field);
    var fe5 = FieldElement.init(5, field);
    var fe6 = FieldElement.init(6, field);
    var fe7 = FieldElement.init(7, field);

    try std.testing.expectEqual(fe0.neg(), fe0);
    try std.testing.expectEqual(fe1.neg(), fe6);
    try std.testing.expectEqual(fe2.neg(), fe5);
    try std.testing.expectEqual(fe3.neg(), fe4);
    try std.testing.expectEqual(fe4.neg(), fe3);
    try std.testing.expectEqual(fe5.neg(), fe2);
    try std.testing.expectEqual(fe6.neg(), fe1);
    try std.testing.expectEqual(fe7.neg(), fe0);
}

test "inverse field elements" {
    var field = Field.init(7);

    var fe0 = field.zero();
    var fe1 = field.one();
    var feminus1 = FieldElement.init(-1, field);
    var fe2 = FieldElement.init(2, field);
    var feminus2 = FieldElement.init(-2, field);
    var fe3 = FieldElement.init(3, field);
    var feminus3 = FieldElement.init(-3, field);
    var fe4 = FieldElement.init(4, field);
    var fe5 = FieldElement.init(5, field);
    var fe6 = FieldElement.init(6, field);
    var fe7 = FieldElement.init(7, field);

    try std.testing.expectEqual(fe0.inv(), fe0);
    try std.testing.expectEqual(fe1.inv(), fe1);
    try std.testing.expectEqual(fe2.inv(), feminus3); 
    try std.testing.expectEqual(fe3.inv(), feminus2); 
    try std.testing.expectEqual(fe4.inv(), fe2);
    try std.testing.expectEqual(fe5.inv(), fe3);
    try std.testing.expectEqual(fe6.inv(), feminus1); 
    try std.testing.expectEqual(fe7.inv(), fe0);
}

test "exponentiating field elements" {
    var field = Field.init(19);

    var fe1 = FieldElement.init(1, field);
    var fe5 = FieldElement.init(5, field);
    var fe6 = FieldElement.init(6, field);

    try std.testing.expectEqual(fe5, fe5.exp(1));
    try std.testing.expectEqual(fe6, fe5.exp(2));
    try std.testing.expectEqual(fe1, fe5.exp(117));
}

pub const Field = struct {
    const Self = @This();
    const GENERATOR_VALUE: FieldElement.TValue = 85408008396924667383611388730472331217;

    const TPrime = i256;

    prime: TPrime,


    pub fn init(prime: TPrime) Self {
        return Self {
            .prime = prime,
        };
    }

    pub fn zero(self: Self) FieldElement {
        return FieldElement.init(0, self);
    }

    pub fn one(self: Self) FieldElement {
        return FieldElement.init(1, self);
    }

    pub fn add(self: Self, left: FieldElement, right: FieldElement) FieldElement {
        return FieldElement.init(@mod((left.value + right.value), self.prime), self);
    }

    pub fn multiply(self: Self, left: FieldElement, right: FieldElement) FieldElement {
        return FieldElement.init(@mod((left.value * right.value), self.prime), self);
    }

    // TODO savil. I don't know why Alan added the prime field here.
    //
    // (5 - 2) % 7 == 3 % 7 == 3
    // (7 + 5 - 2)% 7 == 10 % 7 == 3
    //
    // (2 - 5) % 7 == (-3) % 7 == 4 % 7 = 4
    // (7 + 2 - 5) % 7 == 4 % 7
    //
    // (2 - 20) % 7 == -18 % 7 = -4 % 7 = 3
    // (7 + 2 - 20) % 7 == (-11) % 7 == -4 % 7 = 3
    pub fn subtract(self: Self, left: FieldElement, right: FieldElement) FieldElement {
        return FieldElement.init(@mod((self.prime + left.value - right.value), self.prime), self);
    }

    pub fn negate(self: Self, operand: FieldElement) FieldElement {
        return FieldElement.init(@mod((self.prime - operand.value), self.prime), self);        
    }

    pub fn inverse(self: Self, operand: FieldElement) !FieldElement {
        var result = try xgcd(@intCast(i64, operand.value), @intCast(i64, self.prime));
        return FieldElement.init(result.a, self);
    }

    pub fn divide(self: Self, left: FieldElement, right: FieldElement) !FieldElement {
        var result = try xgcd(@intCast(i64, right.value), @intCast(i64, self.prime));
        return FieldElement.init((left.value * @intCast(u64, result.a)) % self.p, self);
    }

    pub fn generator(self: Self) FieldElement {
        assert(self.prime == 1 + 407 * ( 1 << 119 )); // , "Do not know generator for other fields beyond 1+407*2^119");
        return FieldElement.init(Field.GENERATOR_VALUE, self);
    }

    pub fn primitiveNthRoot(self: Self, n: i64) FieldElement {
        assert(self.prime == 1 + 407 * ( 1 << 119 )); //, "Cannot return root of unity for unknown field");
        assert(n <= (1 << 119)); //  "Field does not have nth root of unity where n > 2^119.");
        assert((n & (n-1)) == 0); // "Field does not have nth root of unity where n is not power of two.");
        var root = FieldElement(85408008396924667383611388730472331217, self);
        var order = 1 << 119;
        while (order != n) {
            root = root^2;
            order = order/2;
        }
        return root;
    }
};


test "generator for field elements" {
    var field = Field.init(1 + 407 * ( 1 << 119));

    try std.testing.expectEqual(Field.GENERATOR_VALUE, field.generator().value);
}

/// xgcd is the extended euclidean algorithm for computing multiplicative inverses:
/// for input a and b, returns their greatest common divisor g, 
/// along with Bezout coefficients a,b such that ax + by = g.
/// https://aszepieniec.github.io/stark-anatomy/basic-tools#finite-field-arithmetic
///
/// The multiplicative inverse of A is A^-1 such that:
/// A * A^-1 = 1 mod C or equivalently (A * A^-1) mod C = 1
/// https://www.khanacademy.org/computing/computer-science/cryptography/modarithmetic/a/modular-inverses
///
const xgcdResult = struct {a: i64, b: i64, g: i64};
fn xgcd(x_in: i64, y_in: i64) !xgcdResult {
    var x = x_in;
    var y = y_in;
    var a: i64 = 1;
    var new_a: i64 = 0;
    var b: i64 = 0;
    var new_b: i64 = 1;

    while (y != 0) {
        var quotient = try std.math.divFloor(@TypeOf(x), x, y);

        var old_x = x;
        x = y;
        y = old_x - quotient * y; // this is the remainder

        var old_a = a;
        a = new_a;
        new_a = old_a - quotient * new_a;

        var old_b = b;
        b = new_b;
        new_b = old_b - quotient * new_b;
    }

    return xgcdResult{.a = a, .b = b, .g = x};
}

test "xgcd" {
    var result = try xgcd(1, 1);
    try std.testing.expectEqual(@as(i64, 1), result.g);

    result = try xgcd(1, 2);
    try std.testing.expectEqual(@as(i64, 1), result.g);

    result = try xgcd(10, 4);
    try std.testing.expectEqual(@as(i64, 2), result.g);

    // example from:
    // https://www.khanacademy.org/computing/computer-science/cryptography/modarithmetic/a/the-euclidean-algorithm
    result = try xgcd(270, 192);
    try std.testing.expectEqual(@as(i64, 6), result.g);
}