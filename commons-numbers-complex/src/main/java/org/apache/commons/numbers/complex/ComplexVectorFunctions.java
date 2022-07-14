package org.apache.commons.numbers.complex;

public class ComplexVectorFunctions {

        /**
         * A complex number representing one.
         *
         * <p>\( (1 + i 0) \).
         */
        public static final ComplexConstructor<Complex> ONE = (x, y) -> Complex.ofCartesian(1, 0);
        /**
         * A complex number representing zero.
         *
         * <p>\( (0 + i 0) \).
         */
        public static final ComplexConstructor<Complex> ZERO = (x, y) -> Complex.ofCartesian(0, 0);

        /** A complex number representing {@code NaN + i NaN}. */
        public static final ComplexConstructor<Complex> NAN = (x, y) -> Complex.ofCartesian(Double.NaN, Double.NaN);
        /** &pi;/2. */
        private static final double PI_OVER_2 = 0.5 * Math.PI;
        /** &pi;/4. */
        private static final double PI_OVER_4 = 0.25 * Math.PI;
        /** Natural logarithm of 2 (ln(2)). */
        private static final double LN_2 = Math.log(2);
        /** Base 10 logarithm of 10 divided by 2 (log10(e)/2). */
        private static final double LOG_10E_O_2 = Math.log10(Math.E) / 2;
        /** Base 10 logarithm of 2 (log10(2)). */
        private static final double LOG10_2 = Math.log10(2);
        /** {@code 1/2}. */
        private static final double HALF = 0.5;
        /** {@code sqrt(2)}. */
        private static final double ROOT2 = 1.4142135623730951;
        /** {@code 1.0 / sqrt(2)}.
         * This is pre-computed to the closest double from the exact result.
         * It is 1 ULP different from 1.0 / Math.sqrt(2) but equal to Math.sqrt(2) / 2.
         */
        private static final double ONE_OVER_ROOT2 = 0.7071067811865476;
        /** The bit representation of {@code -0.0}. */
        private static final long NEGATIVE_ZERO_LONG_BITS = Double.doubleToLongBits(-0.0);
        /** Exponent offset in IEEE754 representation. */
        private static final int EXPONENT_OFFSET = 1023;
        /**
         * Largest double-precision floating-point number such that
         * {@code 1 + EPSILON} is numerically equal to 1. This value is an upper
         * bound on the relative error due to rounding real numbers to double
         * precision floating-point numbers.
         *
         * <p>In IEEE 754 arithmetic, this is 2<sup>-53</sup>.
         * Copied from o.a.c.numbers.Precision.
         *
         * @see <a href="http://en.wikipedia.org/wiki/Machine_epsilon">Machine epsilon</a>
         */
        private static final double EPSILON = Double.longBitsToDouble((EXPONENT_OFFSET - 53L) << 52);
        /** Mask to remove the sign bit from a long. */
        private static final long UNSIGN_MASK = 0x7fff_ffff_ffff_ffffL;
        /** Mask to extract the 52-bit mantissa from a long representation of a double. */
        private static final long MANTISSA_MASK = 0x000f_ffff_ffff_ffffL;
        /** The multiplier used to split the double value into hi and low parts. This must be odd
         * and a value of 2^s + 1 in the range {@code p/2 <= s <= p-1} where p is the number of
         * bits of precision of the floating point number. Here {@code s = 27}.*/
        private static final double MULTIPLIER = 1.34217729E8;

        /**
         * Crossover point to switch computation for asin/acos factor A.
         * This has been updated from the 1.5 value used by Hull et al to 10
         * as used in boost::math::complex.
         * @see <a href="https://svn.boost.org/trac/boost/ticket/7290">Boost ticket 7290</a>
         */
        private static final double A_CROSSOVER = 10.0;
        /** Crossover point to switch computation for asin/acos factor B. */
        private static final double B_CROSSOVER = 0.6471;
        /**
         * The safe maximum double value {@code x} to avoid loss of precision in asin/acos.
         * Equal to sqrt(M) / 8 in Hull, et al (1997) with M the largest normalised floating-point value.
         */
        private static final double SAFE_MAX = Math.sqrt(Double.MAX_VALUE) / 8;
        /**
         * The safe minimum double value {@code x} to avoid loss of precision/underflow in asin/acos.
         * Equal to sqrt(u) * 4 in Hull, et al (1997) with u the smallest normalised floating-point value.
         */
        private static final double SAFE_MIN = Math.sqrt(Double.MIN_NORMAL) * 4;
        /**
         * The safe maximum double value {@code x} to avoid loss of precision in atanh.
         * Equal to sqrt(M) / 2 with M the largest normalised floating-point value.
         */
        private static final double SAFE_UPPER = Math.sqrt(Double.MAX_VALUE) / 2;
        /**
         * The safe minimum double value {@code x} to avoid loss of precision/underflow in atanh.
         * Equal to sqrt(u) * 2 with u the smallest normalised floating-point value.
         */
        private static final double SAFE_LOWER = Math.sqrt(Double.MIN_NORMAL) * 2;
        /** The safe maximum double value {@code x} to avoid overflow in sqrt. */
        private static final double SQRT_SAFE_UPPER = Double.MAX_VALUE / 8;
        /**
         * A safe maximum double value {@code m} where {@code e^m} is not infinite.
         * This can be used when functions require approximations of sinh(x) or cosh(x)
         * when x is large using exp(x):
         * <pre>
         * sinh(x) = (e^x - e^-x) / 2 = sign(x) * e^|x| / 2
         * cosh(x) = (e^x + e^-x) / 2 = e^|x| / 2 </pre>
         *
         * <p>This value can be used to approximate e^x using a product:
         *
         * <pre>
         * e^x = product_n (e^m) * e^(x-nm)
         * n = (int) x/m
         * e.g. e^2000 = e^m * e^m * e^(2000 - 2m) </pre>
         *
         * <p>The value should be below ln(max_value) ~ 709.783.
         * The value m is set to an integer for less error when subtracting m and chosen as
         * even (m=708) as it is used as a threshold in tanh with m/2.
         *
         * <p>The value is used to compute e^x multiplied by a small number avoiding
         * overflow (sinh/cosh) or a small number divided by e^x without underflow due to
         * infinite e^x (tanh). The following conditions are used:
         * <pre>
         * 0.5 * e^m * Double.MIN_VALUE * e^m * e^m = Infinity
         * 2.0 / e^m / e^m = 0.0 </pre>
         */
        private static final double SAFE_EXP = 708;
        /**
         * The value of Math.exp(SAFE_EXP): e^708.
         * To be used in overflow/underflow safe products of e^m to approximate e^x where x > m.
         */
        private static final double EXP_M = Math.exp(SAFE_EXP);

        /** 54 shifted 20-bits to align with the exponent of the upper 32-bits of a double. */
        private static final int EXP_54 = 0x36_00000;
        /** Represents an exponent of 500 in unbiased form shifted 20-bits to align with the upper 32-bits of a double. */
        private static final int EXP_500 = 0x5f3_00000;
        /** Represents an exponent of 1024 in unbiased form (infinite or nan)
         * shifted 20-bits to align with the upper 32-bits of a double. */
        private static final int EXP_1024 = 0x7ff_00000;
        /** Represents an exponent of -500 in unbiased form shifted 20-bits to align with the upper 32-bits of a double. */
        private static final int EXP_NEG_500 = 0x20b_00000;
        /** 2^600. */
        private static final double TWO_POW_600 = 0x1.0p+600;
        /** 2^-600. */
        private static final double TWO_POW_NEG_600 = 0x1.0p-600;

        /** Serializable version identifier. */
        private static final long serialVersionUID = 20180201L;

        /**
         * The size of the buffer for {@link #toString()}.
         *
         * <p>The longest double will require a sign, a maximum of 17 digits, the decimal place
         * and the exponent, e.g. for max value this is 24 chars: -1.7976931348623157e+308.
         * Set the buffer size to twice this and round up to a power of 2 thus
         * allowing for formatting characters. The size is 64.
         */
        private static final int TO_STRING_SIZE = 64;
        /** The minimum number of characters in the format. This is 5, e.g. {@code "(0,0)"}. */
        private static final int FORMAT_MIN_LEN = 5;
        /** {@link #toString() String representation}. */
        private static final char FORMAT_START = '(';
        /** {@link #toString() String representation}. */
        private static final char FORMAT_END = ')';
        /** {@link #toString() String representation}. */
        private static final char FORMAT_SEP = ',';
        /** The minimum number of characters before the separator. This is 2, e.g. {@code "(0"}. */
        private static final int BEFORE_SEP = 2;


        /**
         * Check that a value is negative. It must meet all the following conditions:
         * <ul>
         *  <li>it is not {@code NaN},</li>
         *  <li>it is negative signed,</li>
         * </ul>
         *
         * <p>Note: This is true for negative zero.</p>
         *
         * @param d Value.
         * @return {@code true} if {@code d} is negative.
         */
        static boolean negative(double d) {
            return d < 0 || Double.doubleToLongBits(d) == NEGATIVE_ZERO_LONG_BITS;
        }

        /**
         * Change the sign of the magnitude based on the signed value.
         *
         * <p>If the signed value is negative then the result is {@code -magnitude}; otherwise
         * return {@code magnitude}.
         *
         * <p>A signed value of {@code -0.0} is treated as negative. A signed value of {@code NaN}
         * is treated as positive.
         *
         * <p>This is not the same as {@link Math#copySign(double, double)} as this method
         * will change the sign based on the signed value rather than copy the sign.
         *
         * @param magnitude the magnitude
         * @param signedValue the signed value
         * @return magnitude or -magnitude.
         * @see #negative(double)
         */
        private static double changeSign(double magnitude, double signedValue) {
            return negative(signedValue) ? -magnitude : magnitude;
        }

        /**
         * Returns a scale suitable for use with {@link Math#scalb(double, int)} to normalise
         * the number to the interval {@code [1, 2)}.
         *
         * <p>The scale is typically the largest unbiased exponent used in the representation of the
         * two numbers. In contrast to {@link Math#getExponent(double)} this handles
         * sub-normal numbers by computing the number of leading zeros in the mantissa
         * and shifting the unbiased exponent. The result is that for all finite, non-zero,
         * numbers {@code a, b}, the magnitude of {@code scalb(x, -getScale(a, b))} is
         * always in the range {@code [1, 2)}, where {@code x = max(|a|, |b|)}.
         *
         * <p>This method is a functional equivalent of the c function ilogb(double) adapted for
         * two input arguments.
         *
         * <p>The result is to be used to scale a complex number using {@link Math#scalb(double, int)}.
         * Hence the special case of both zero arguments is handled using the return value for NaN
         * as zero cannot be scaled. This is different from {@link Math#getExponent(double)}
         * or {@link #getMaxExponent(double, double)}.
         *
         * <p>Special cases:
         *
         * <ul>
         * <li>If either argument is NaN or infinite, then the result is
         * {@link Double#MAX_EXPONENT} + 1.
         * <li>If both arguments are zero, then the result is
         * {@link Double#MAX_EXPONENT} + 1.
         * </ul>
         *
         * @param a the first value
         * @param b the second value
         * @return The maximum unbiased exponent of the values to be used for scaling
         * @see Math#getExponent(double)
         * @see Math#scalb(double, int)
         * @see <a href="http://www.cplusplus.com/reference/cmath/ilogb/">ilogb</a>
         */
        private static int getScale(double a, double b) {
            // Only interested in the exponent and mantissa so remove the sign bit
            final long x = Double.doubleToRawLongBits(a) & UNSIGN_MASK;
            final long y = Double.doubleToRawLongBits(b) & UNSIGN_MASK;
            // Only interested in the maximum
            final long bits = Math.max(x, y);
            // Get the unbiased exponent
            int exp = ((int) (bits >>> 52)) - EXPONENT_OFFSET;

            // No case to distinguish nan/inf
            // Handle sub-normal numbers
            if (exp == Double.MIN_EXPONENT - 1) {
                // Special case for zero, return as nan/inf to indicate scaling is not possible
                if (bits == 0) {
                    return Double.MAX_EXPONENT + 1;
                }
                // A sub-normal number has an exponent below -1022. The amount below
                // is defined by the number of shifts of the most significant bit in
                // the mantissa that is required to get a 1 at position 53 (i.e. as
                // if it were a normal number with assumed leading bit)
                final long mantissa = bits & MANTISSA_MASK;
                exp -= Long.numberOfLeadingZeros(mantissa << 12);
            }
            return exp;
        }

        /**
         * Returns the largest unbiased exponent used in the representation of the
         * two numbers. Special cases:
         *
         * <ul>
         * <li>If either argument is NaN or infinite, then the result is
         * {@link Double#MAX_EXPONENT} + 1.
         * <li>If both arguments are zero or subnormal, then the result is
         * {@link Double#MIN_EXPONENT} -1.
         * </ul>
         *
         * <p>This is used by {@link #} as
         * a simple detection that a number may overflow if multiplied
         * by a value in the interval [1, 2).
         *
         * @param a the first value
         * @param b the second value
         * @return The maximum unbiased exponent of the values.
         * @see Math#getExponent(double)
         * @see #(double, double, double, double)
         */
        static int getMaxExponent(double a, double b) {
            // This could return:
            // Math.getExponent(Math.max(Math.abs(a), Math.abs(b)))
            // A speed test is required to determine performance.
            return Math.max(Math.getExponent(a), Math.getExponent(b));
        }

        /**
         * Checks if both x and y are in the region defined by the minimum and maximum.
         *
         * @param x x value.
         * @param y y value.
         * @param min the minimum (exclusive).
         * @param max the maximum (exclusive).
         * @return true if inside the region.
         */
        private static boolean inRegion(double x, double y, double min, double max) {
            return (x < max) && (x > min) && (y < max) && (y > min);
        }

        /**
         * Check that a value is positive infinity. Used to replace {@link Double#isInfinite()}
         * when the input value is known to be positive (i.e. in the case where it has been
         * set using {@link Math#abs(double)}).
         *
         * @param d Value.
         * @return {@code true} if {@code d} is +inf.
         */
        private static boolean isPosInfinite(double d) {
            return d == Double.POSITIVE_INFINITY;
        }

        public static ComplexVector conj(ComplexVector input, ComplexVector out) {
            final int len = input.size();
            for (int i = 0; i < len; i++) {
                ComplexDouble c = input.get(i);
                out = arrayConj(c.real(), c.imag(), (ComplexVectorResult<ComplexVector>) out, i);
            }
            return out;
        }

        public static ComplexVector arrayConj(double r, double i,
                                              ComplexVectorResult<ComplexVector> resultConsumer, int index) {
            return resultConsumer.setValue(index, r, -i);
        }

}
