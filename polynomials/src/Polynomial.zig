// TODO savil.
// - interpolateDomain: write tests
// - zerofierDomain: implement
// - scale: implement
// - testColinearity: implement




const std = @import("std");
const testing = std.testing;
const Allocator = std.mem.Allocator;
// TODO savil. Should this import the library instead of the file?
// TODO savil. Rename FieldElement.zig to FiniteFields.zig?
const FiniteFields = @import("finite-fields");
const FieldElement = FiniteFields.FieldElement;
const Field = FiniteFields.Field;




// move to an integer library
/// Returns a slice of integers from start to end (inclusive) that can be iterated over
/// from: https://gist.github.com/hazeycode/e7e2d81ea2b5b9137502cfe04541080e
///
/// TODO savil: how do slices work in memory? In regular code if I return an array
/// residing inside a function then that memory would get cleaned up upon function return.
pub fn range(comptime start: comptime_int, comptime end: comptime_int) []comptime_int {
    const d: isize = end - start;
    const d_norm = if (d < 0) -1 else if (d > 0) 1 else 0;
    const len = (try std.math.absInt(d)) + 1;
    comptime var i = 0;
    comptime var n = start;
    comptime var res: [len]comptime_int = .{undefined} ** len;
    inline while (i < len) : (i += 1) {
        res[i] = n;
        n += d_norm;
    }
    return res[0..len];
}

/// Returns a slice of integers from start (inclusive) to end (exclusive)
pub fn rangeEx(comptime start: comptime_int, comptime end: comptime_int) []comptime_int {
    return range(start, end-1);
}

// TODO savil. If ziglang had varargs, then we could introduce this syntax for constructing
// the coefficients of the polynomial which are an array of FieldElements.
//
// Right now, for a polynomial like 3x^2+4x+5, the elemenets in
// `[]FieldElements{fe5, fe4, fe3}` are in reverse order. 
// I'd much prefer the ergonomics of `coefficients(fe3, fe4, fe5)`.
//
// fn coefficients(...elems: FieldElement) []FieldElement {
// }


const Polynomial = struct {
    const Self = @This();

    // we need an allocator to store the coefficients in the heap
    allocator: *const Allocator,
    coefficients: []FieldElement,
    
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

        var idx: usize = 0;
        while (idx <= self_degree.?) {
            const self_coeff = self.coefficients[idx];
            if (!self_coeff.eq(other.coefficients[idx])) {
                return false;
            }
            idx += 1;
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

    // Almost certainly a faster way to do this.
    fn exp(self: *Self, exponent: i64) !*Self { 
        if (self.isZero()) {
            var buf = try self.allocator.alloc(FieldElement, 0);
            self.allocator.free(self.coefficients);
            self.coefficients = buf;
            return self;
        }

        var polyOne = try Polynomial.init(self.allocator, &[_]FieldElement{self.coefficients[0].field.one()});
        defer polyOne.deinit();
        var result = &polyOne;

        var powOfTwo: i64 = 1;
        var selfCopy = try Polynomial.init(self.allocator, self.coefficients); 
        defer selfCopy.deinit();
        var valuePowOfTwo = &selfCopy;

        while (powOfTwo <= exponent) { 

            if (powOfTwo & exponent == powOfTwo) {
                result = try result.mul(valuePowOfTwo); 
            }

            valuePowOfTwo = try valuePowOfTwo.mul(valuePowOfTwo);

            powOfTwo = powOfTwo << 1;
        }

        // Copy over the result coefficients into self. Can probably eliminate this copy
        // by directly using self, instead of result.
        try self.writeCoeffs(result.coefficients);

        return self;
    }

    const divResult = struct {
        quotient: Polynomial,
        remainder: Polynomial,

        fn init(q: Polynomial, r: Polynomial) divResult {
            return divResult {
                .quotient = q,
                .remainder = r,
            };
        }

        // for convenience in test cleanup
        fn deinit(self: *divResult) void {
            self.quotient.deinit();
            self.remainder.deinit();
        }

        fn eq(self: *divResult, other: *const divResult) bool {
            return self.quotient.eq(&other.quotient) and 
                self.remainder.eq(&other.remainder);
        }
    };

    /// NOTE: we place div as a static helper method
    //
    // This generally follows the division algorithm from school:
    // keep dividing the dividend by the divisor, until the remainder's highest degree 
    // is less than that of the divisor.
    //
    // ownership of the quotient and remainder is passed on to the caller.
    fn div(allocator: *const Allocator, numerator: *Polynomial, denominator: *Polynomial) !divResult {
        var denom_degree = denominator.degree();
        if (denom_degree == null) {
            return error.DivisionByZero;
        }
        var num_degree = numerator.degree();
        if (num_degree == null or num_degree.? < denom_degree.?) {
            var empty_poly = try Polynomial.init(allocator, &[_]FieldElement{});
            // copy the numerator so that the caller owns the new value
            return divResult.init(empty_poly, try numerator.copy());
        }

        var field = denominator.coefficients[0].field;
        // remainder:
        var remainder = try Polynomial.init(allocator, numerator.coefficients);

        var quotient_coeffs_size = num_degree.? - denom_degree.? + 1;
        var quotient_coeffs = try allocator.alloc(FieldElement, quotient_coeffs_size);
        defer allocator.free(quotient_coeffs);
        std.mem.set(FieldElement, quotient_coeffs, field.zero());

        var idx: usize = 0;
        while (idx < quotient_coeffs_size): (idx += 1) {
            const rem_degree = remainder.degree();
            if (rem_degree == null or rem_degree.? < denom_degree.?) {
                // The remainder's highest degree is less than the divisor's 
                // highest degree, so we should stop here.
                break;
            }
            
            var coeff = try remainder.leadingCoefficient().?.div(denominator.leadingCoefficient().?);
            var shift = rem_degree.? - denom_degree.?;
            
            // Initialize a subtractee polynomial whose leading coefficient
            // is `coeff` and is of length `shift + 1`
            var subtractee_coeffs = try allocator.alloc(FieldElement, shift + 1);
            defer allocator.free(subtractee_coeffs);
            std.mem.set(FieldElement, subtractee_coeffs, field.zero());
            subtractee_coeffs[subtractee_coeffs.len - 1] = coeff;
            var subtractee = &(try Polynomial.init(allocator, subtractee_coeffs));
            defer subtractee.deinit();
            // Finally, multiply the subtractee by the denominator (i.e. divisor)
            subtractee = try subtractee.mul(denominator);
            
            quotient_coeffs[shift] = coeff;

            _ = try remainder.sub(subtractee);
        }

        var q = try Polynomial.init(allocator, quotient_coeffs);
        return divResult.init(q, remainder);
    }

    fn quot(self: *Self, other: *Polynomial) !*Self {
        var result = try Polynomial.div(self.allocator, self, other);
        defer result.deinit();
        try writeCoeffs(self, result.quotient.coefficients);
        return self;
    }

    fn rem(self: *Self, other: *Polynomial) !*Self {
        var result = try Polynomial.div(self.allocator, self, other);
        defer result.deinit();
        try writeCoeffs(self, result.remainder.coefficients);
        return self;
    }

    fn evaluate(self: *Self, point: FieldElement) FieldElement {

        var x_i = point.field.one();
        var value = point.field.zero();

        for (self.coefficients) |c| {
            value = value.add(c.mul(x_i));
            x_i = x_i.mul(point);
        }
        return value;
    }

    // the caller must free the result value
    //
    // We explicitly pass in an allocator so that the caller is in charge of the memory.
    // Otherwise, it would be weird that a value allocator from an allcoator internal to Polynomial
    // is expected to be freed by the caller using that internal allocator. 
    fn evaluateDomain(self: *Self, a: *const Allocator, domain: []FieldElement) ![]FieldElement {
        var results = try a.alloc(FieldElement, domain.len); 
        var i: usize = 0;
        for (domain) | val | {
            results[i] = self.evaluate(val);
            i += 1;
        }
        return results;
    }

    fn interpolateDomain(
        allocator: *Allocator, 
        domain: []FieldElement, 
        values: []FieldElement,
    ) !Polynomial {
        if (domain.len != values.len) {
            return error.DomainLenIsNotSameAsValuesLen;
        }
        if (domain.len == 0) {
            return error.DomainLenIsZero;
        }

        var field = domain[0].field();
        var acc = try Polynomial.init(allocator, &[_]FieldElement{});
        defer acc.deinit();

        for (rangeEx(0, domain.len)) | i | {
            var product = Polynomial.init(allocator, &[_]FieldElement{values[i]});
            for (rangeEx(0, domain.len)) | j | {
                if (j == 1) {
                    continue;
                }
                // CODE SMELL: we add x inside this inner loop instead of defining it once outside
                // because the Polynomial sub function mutates the polynomial it is called upon.
                // An alternate design may be to have the sub function return a newly allocated polynomial.
                var x = try Polynomial.init(allocator, &[_]FieldElement{field.zero(), field.one()});

                var polyDomainOfJ = Polynomial.init(allocator, &[_]FieldElement{domain[j]});
                var polyXMinusDomainOfJ = x.sub(polyDomainOfJ);

                var inverseOfIMinusJ = (domain[i] - domain[j]).inverse();
                var polyDomainOfIMinusJ = Polynomial.init(allocator, &[_]FieldElement{inverseOfIMinusJ});

                product = product.mul(polyXMinusDomainOfJ.mul(polyDomainOfIMinusJ));

                x.deinit();
                polyDomainOfJ.deinit();
                polyDomainOfIMinusJ.deinit();
            }
            acc = acc.add(product);
            product.deinit();
        }

    }

    fn copy(self: *Self) !Polynomial {
        return try Polynomial.init(self.allocator, self.coefficients);
    }

    fn writeCoeffs(self: *Self, coeffs: []FieldElement) !void {
        self.allocator.free(self.coefficients);
        self.coefficients = try self.allocator.alloc(FieldElement, coeffs.len);
        std.mem.copy(FieldElement, self.coefficients, coeffs);
    }
};

test "interpolateDomain" {
    // const field = Field.init(19);
    // const fe0 = field.zero();
    // const fe1 = field.one();
    // const fe2 = FieldElement.init(2, field);
    // const fe4 = FieldElement.init(4, field);
    // const fe7 = FieldElement.init(7, field);
    // const fe16 = FieldElement.init(16, field);

    // Polynomial.interpolateDomain(&[_]FieldElement{});
}

test "evaluateDomain" {

    const field = Field.init(19);
    //const fe0 = field.zero();
    const fe1 = field.one();
    const fe2 = FieldElement.init(2, field);
    const fe4 = FieldElement.init(4, field);
    const fe7 = FieldElement.init(7, field);
    const fe16 = FieldElement.init(16, field);

    var polyDegTwo = &(try Polynomial.init(&testing.allocator, &[_]FieldElement{fe1, fe2, fe4}));
    defer polyDegTwo.deinit();

    var values = try polyDegTwo.evaluateDomain(&testing.allocator, &[_]FieldElement{fe1, fe4, fe4});
    defer testing.allocator.free(values);

    var expectedResults = [_]FieldElement{fe7, fe16, fe16};
    var i: usize = 0;
    for (values) | value | {
        try testing.expect(value.eq(expectedResults[i]));
        i += 1;
    }
}

test "evaluate" {

    const field = Field.init(19);
    const fe0 = field.zero();
    const fe1 = field.one();
    const fe2 = FieldElement.init(2, field);
    const fe4 = FieldElement.init(4, field);
    const fe9 = FieldElement.init(9, field);
    const fe16 = FieldElement.init(16, field);

    var polyEmpty = &(try Polynomial.init(&testing.allocator, &[_]FieldElement{}));
    defer polyEmpty.deinit();
    
    const resultEmpty = polyEmpty.evaluate(fe4);
    try testing.expect(resultEmpty.eq(fe0));

    // polynomial of degree zero
    var polyDegZero = &(try Polynomial.init(&testing.allocator, &[_]FieldElement{fe2}));
    defer polyDegZero.deinit();

    const resultDegZero = polyDegZero.evaluate(fe4);
    try testing.expect(resultDegZero.eq(fe2));

    // polynomial of degree one
    var polyDegOne = &(try Polynomial.init(&testing.allocator, &[_]FieldElement{fe1, fe2}));
    defer polyDegOne.deinit();

    const resultDegOne = polyDegOne.evaluate(fe4);
    try testing.expect(resultDegOne.eq(fe9));

    // polynomial of degree two
    var polyDegTwo = &(try Polynomial.init(&testing.allocator, &[_]FieldElement{fe1, fe2, fe4}));
    defer polyDegTwo.deinit();

    const resultDegTwo = polyDegTwo.evaluate(fe4);
    // = 1 + 2x + 4x^2 = 1 + (2)(4) + (4)(4^2) =  73 mod 19 = 16
    try testing.expect(resultDegTwo.eq(fe16));
}

test "quotient" {
    const field = Field.init(19);
    const fe0 = field.zero();
    const fe1 = field.one();
    const fe2 = FieldElement.init(2, field);
    const fe4 = FieldElement.init(4, field);

    // Test Case: 2x^2 + 4x + 4 / x + 2 = (quotient = 2x, remainder = 4)

    var dividend = &(try Polynomial.init(&testing.allocator, &[_]FieldElement{fe4, fe4, fe2}));
    defer dividend.deinit();    

    var divisor = &(try Polynomial.init(&testing.allocator, &[_]FieldElement{fe2, fe1}));
    defer divisor.deinit();    

    _ = try dividend.quot(divisor);

    var expectedQuotient = &(try Polynomial.init(&testing.allocator, &[_]FieldElement{fe0, fe2}));
    defer expectedQuotient.deinit();
    try testing.expect(dividend.eq(expectedQuotient));
}

test "remainder" {
    const field = Field.init(19);
    //const fe0 = field.zero();
    const fe1 = field.one();
    const fe2 = FieldElement.init(2, field);
    const fe4 = FieldElement.init(4, field);

    // Test Case: 2x^2 + 4x + 4 / x + 2 = (quotient = 2x, remainder = 4)

    var dividend = &(try Polynomial.init(&testing.allocator, &[_]FieldElement{fe4, fe4, fe2}));
    defer dividend.deinit();    

    var divisor = &(try Polynomial.init(&testing.allocator, &[_]FieldElement{fe2, fe1}));
    defer divisor.deinit();    

    _ = try dividend.rem(divisor);

    var expected_rem = &(try Polynomial.init(&testing.allocator, &[_]FieldElement{fe4}));
    defer expected_rem.deinit();
    try testing.expect(dividend.eq(expected_rem));
}

test "div-empty" {
    const emptyPoly = &(try Polynomial.init(&testing.allocator, &[_]FieldElement{}));
    defer emptyPoly.deinit();

    var emptyPoly2 = &(try emptyPoly.copy());
    defer emptyPoly2.deinit();

    // Test Case: 0 / 0 = null
    var result = Polynomial.div(&testing.allocator, emptyPoly, emptyPoly2);
    try testing.expectError(error.DivisionByZero, result);
}

test "div-zero" {

    const field = Field.init(19);
    const fe1 = field.one();
    const fe2 = FieldElement.init(2, field);
    
    const emptyPoly = &(try Polynomial.init(&testing.allocator, &[_]FieldElement{}));
    defer emptyPoly.deinit();

    // Test Case: 0 / x + 2 = (quotient = 0, remainder = 0)
    try testDivHelper(
        &[_]FieldElement{}, // dividend
        &[_]FieldElement{fe2, fe1}, // divisor
        &[_]FieldElement{}, // quotient
        &[_]FieldElement{}, // remainder
    );
}

test "div-by-self" {
    const field = Field.init(19);
    const fe1 = field.one();

    // Test case: x + 1 / x + 1 = (quotient = 0, remainder = 1)
    try testDivHelper(
        &[_]FieldElement{fe1, fe1}, // dividend
        &[_]FieldElement{fe1, fe1}, // divisor
        &[_]FieldElement{fe1}, // quotient
        &[_]FieldElement{}, // remainder
    );
}

test "div-a-bit-complicated" {
    const field = Field.init(19);
    const fe0 = field.zero();
    const fe1 = field.one();
    const fe2 = FieldElement.init(2, field);
    const fe4 = FieldElement.init(4, field);

    // Test Case: 2x^2 + 4x + 4 / x + 2 = (quotient = 2x, remainder = 4)
    try testDivHelper(
        &[_]FieldElement{fe4, fe4, fe2}, // dividend
        &[_]FieldElement{fe2, fe1}, // divisor
        &[_]FieldElement{fe0, fe2}, // quotient
        &[_]FieldElement{fe4}, // remainder
    );
}

fn testDivHelper(
    dividendCoeffs: []FieldElement, 
    divisorCoeffs: []FieldElement, 
    quotientCoeffs: []FieldElement,
    remainderCoeffs: []FieldElement,
) !void {
    const dividend = &(try Polynomial.init(&testing.allocator, dividendCoeffs));
    defer dividend.deinit();

    const divisor = &(try Polynomial.init(&testing.allocator, divisorCoeffs));
    defer divisor.deinit();

    const quotient = (try Polynomial.init(&testing.allocator, quotientCoeffs));
    defer quotient.deinit();

    const remainder = (try Polynomial.init(&testing.allocator, remainderCoeffs));
    defer remainder.deinit();

    var result = try Polynomial.div(&testing.allocator, dividend, divisor);
    defer result.deinit();

    var expected = Polynomial.divResult{.quotient = quotient, .remainder = remainder};
    try testing.expect(result.eq(&expected));
}

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

test "degrees-trailing-zeros" {
    const field = Field.init(19);
    const fe0 = field.zero();
    const fe2 = FieldElement.init(2, field);
    const polynomial = try Polynomial.init(&testing.allocator, &[_]FieldElement{fe0, fe2, fe0, fe0});
    defer polynomial.deinit();
    try testing.expect(polynomial.degree().? == 1);
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

test "exponential" {
    const field = Field.init(19);
    const fe1 = FieldElement.init(1, field);
    const fe2 = FieldElement.init(2, field);
    const fe3 = FieldElement.init(3, field);
    const fe4 = FieldElement.init(4, field);
    const fe5 = FieldElement.init(5, field);
    const fe6 = FieldElement.init(6, field);
    const fe8 = FieldElement.init(8, field);
    const fe11 = FieldElement.init(11, field);
    const fe17 = FieldElement.init(17, field);

    var polynomial = try Polynomial.init(&testing.allocator, &[_]FieldElement{fe2, fe5});
    defer polynomial.deinit();

    var test_poly1 = try Polynomial.init(&testing.allocator, &[_]FieldElement{fe2, fe5});
    defer test_poly1.deinit();

    var exp0_poly = try test_poly1.exp(0);
    // sanity test that the identical_polynomial was also changed.
    try testing.expect(exp0_poly.eq(&test_poly1));
    // not equal because polynomial is now field.one()
    try testing.expect(exp0_poly.neq(&polynomial));

    // Raise to the power of 1
    try testExpHelper(&[_]FieldElement{fe2, fe5}, 1);    

    // Raise to power of 2
    try testExpHelper(&[_]FieldElement{fe4, fe1, fe6}, 2);

    // Raise to power of 3 
    //(5x+2)^3 = (25x^2 + 20x + 4)(5x+2) mod 19 = 125x^3+100x^2+20x+50x^2+40x+8 = 11x^3+5x^2+1x+12x^2+2x+8=11x^3+17x^2+3x+8
    try testExpHelper(&[_]FieldElement{fe8, fe3, fe17, fe11}, 3);

    // Raie to the power of 4   
    // (11x^3+17x^2+3x+8)(5x+2) = feeling sleepy
    //try testExpHelper(&[_]FieldElement{fe8, fe3, fe17, fe11}, 4);
}

fn testExpHelper(expected_coefficients: []FieldElement, exp: i64) !void {
    const field = Field.init(19);
    const fe2 = FieldElement.init(2, field);
    const fe5 = FieldElement.init(5, field);

    var test_poly3 = try Polynomial.init(&testing.allocator, &[_]FieldElement{fe2, fe5});
    defer test_poly3.deinit();
    _ = try test_poly3.exp(exp);
    
    var expected_poly3 = try Polynomial.init(&testing.allocator, expected_coefficients);
    defer expected_poly3.deinit();
    try testing.expect(expected_poly3.eq(&test_poly3));
}