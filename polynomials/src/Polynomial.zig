const std = @import("std");
const testing = std.testing;
// TODO savil. Should this import the library instead of the file?
// TODO savil. Rename FieldElement.zig to FiniteFields.zig?
const FiniteFields = @import("finite-fields");
const FieldElement = FiniteFields.FieldElement;
const Field = FiniteFields.Field;


const Polynomial = struct {
    const Self = @This();

    coefficients: []const FieldElement,

    pub fn init(coefficients: []const FieldElement) Self {
        return Self {
            .coefficients = coefficients,
        };
    }
};

test "init polynomial" {
    const field = Field.init(19);
    const fe = FieldElement.init(2, field);
    const fes = [_]FieldElement{fe};
    const polynomial = Polynomial.init(fes[0..]);

    _ = polynomial;
    try testing.expect(true); // just want the above to compile.
}

