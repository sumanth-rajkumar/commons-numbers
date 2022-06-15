/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package org.apache.commons.numbers.complex;


import java.util.function.DoubleBinaryOperator;
import java.util.function.DoubleUnaryOperator;

public final class ComplexFunctions {

    /**
     * A complex number representing one.
     *
     * <p>\( (1 + i 0) \).
     */
    public static final ComplexResult<Complex> ONE = (x, y) -> Complex.ofCartesian(1, 0);
    /**
     * A complex number representing zero.
     *
     * <p>\( (0 + i 0) \).
     */
    public static final ComplexResult<Complex> ZERO = (x, y) -> Complex.ofCartesian(0, 0);

    /** A complex number representing {@code NaN + i NaN}. */
    public static final ComplexResult<Complex> NAN = (x, y) -> Complex.ofCartesian(Double.NaN, Double.NaN);
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

    private ComplexFunctions() {

    }
    public static <R> R negate(double r, double i, ComplexResult<R> result) {
        return result.apply(-r, -i);
    }
    public static <R> R multiplyIma(double r, double i, ComplexResult<R> result) {
        return result.apply(-1 * i, r);
    }
    /**
     * Returns {@code true} if either real or imaginary component of the complex number is infinite.
     *
     * <p>Note: A complex number with at least one infinite part is regarded
     * as an infinity (even if its other part is a NaN).
     *
     * @return {@code true} if this instance contains an infinite value.
     * @see Double#isInfinite(double)
     */
    public static boolean isInfinite(double real, double imaginary) {
        return Double.isInfinite(real) || Double.isInfinite(imaginary);
    }
    public static <R> R proj(double r, double i, ComplexResult<R> result) {
        if (isInfinite(r, i)) {
            return result.apply(Double.POSITIVE_INFINITY, Math.copySign(0.0, i));
        }
        return result.apply(r, i);
    }
    public static double arg(double r, double i) {
        // Delegate
        return Math.atan2(r, i);
    }
    /**
     * Returns the squared norm value of this complex number. This is also called the absolute
     * square.
     *
     * <p>\[ \text{norm}(x + i y) = x^2 + y^2 \]
     *
     * <p>If either component is infinite then the result is positive infinity. If either
     * component is NaN and this is not {@link #isInfinite(double, double) infinite} then the result is NaN.
     *
     * <p>Note: This method may not return the same value as the square of {@link #/abs()} as
     * that method uses an extended precision computation.
     *
     * <p>{@code norm()} can be used as a faster alternative than {@code abs()} for ranking by
     * magnitude. If used for ranking any overflow to infinity will create an equal ranking for
     * values that may be still distinguished by {@code abs()}.
     *
     * @return The square norm value.
     * @see #isInfinite(double, double)
     * @see #//isNaN()
     * @see #//abs()
     * @see <a href="http://mathworld.wolfram.com/AbsoluteSquare.html">Absolute square</a>
     */
    public static double norm(double real, double imaginary) {
        if (isInfinite(real, imaginary)) {
            return Double.POSITIVE_INFINITY;
        }
        return real * real + imaginary * imaginary;
    }



    /**
     * Returns the absolute value of the complex number.
     * <pre>abs(x + i y) = sqrt(x^2 + y^2)</pre>
     *
     * <p>This should satisfy the special cases of the hypot function in ISO C99 F.9.4.3:
     * "The hypot functions compute the square root of the sum of the squares of x and y,
     * without undue overflow or underflow."
     *
     * <ul>
     * <li>hypot(x, y), hypot(y, x), and hypot(x, −y) are equivalent.
     * <li>hypot(x, ±0) is equivalent to |x|.
     * <li>hypot(±∞, y) returns +∞, even if y is a NaN.
     * </ul>
     *
     * <p>This method is called by all methods that require the absolute value of the complex
     * number, e.g. abs(), sqrt() and log().
     *
     * @param real Real part.
     * @param imaginary Imaginary part.
     * @return The absolute value.
     */
    public static double abs(double real, double imaginary) {
        // Specialised implementation of hypot.
        // See NUMBERS-143
        return hypot(real, imaginary);
    }
    /**
     * Returns {@code sqrt(x^2 + y^2)} without intermediate overflow or underflow.
     *
     * <p>Special cases:
     * <ul>
     * <li>If either argument is infinite, then the result is positive infinity.
     * <li>If either argument is NaN and neither argument is infinite, then the result is NaN.
     * </ul>
     *
     * <p>The computed result is expected to be within 1 ulp of the exact result.
     *
     * <p>This method is a replacement for {@link Math#hypot(double, double)}. There
     * will be differences between this method and {@code Math.hypot(double, double)} due
     * to the use of a different algorithm to compute the high precision sum of
     * {@code x^2 + y^2}. This method has been tested to have a lower maximum error from
     * the exact result; any differences are expected to be 1 ULP indicating a rounding
     * change in the sum.
     *
     * <p>JDK9 ported the hypot function to Java for bug JDK-7130085 due to the slow performance
     * of the method as a native function. Benchmarks of the Complex class for functions that
     * use hypot confirm this is slow pre-Java 9. This implementation outperforms the new faster
     * {@code Math.hypot(double, double)} on JDK 11 (LTS). See the Commons numbers examples JMH
     * module for benchmarks. Comparisons with alternative implementations indicate
     * performance gains are related to edge case handling and elimination of an unpredictable
     * branch in the computation of {@code x^2 + y^2}.
     *
     * <p>This port was adapted from the "Freely Distributable Math Library" hypot function.
     * This method only (and not invoked methods within) is distributed under the terms of the
     * original notice as shown below:
     * <pre>
     * ====================================================
     * Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
     *
     * Developed at SunSoft, a Sun Microsystems, Inc. business.
     * Permission to use, copy, modify, and distribute this
     * software is freely granted, provided that this notice
     * is preserved.
     * ====================================================
     * </pre>
     *
     * <p>Note: The fdlibm c code makes use of the language ability to read and write directly
     * to the upper and lower 32-bits of the 64-double. The function performs
     * checking on the upper 32-bits for the magnitude of the two numbers by accessing
     * the exponent and 20 most significant bits of the mantissa. These upper bits
     * are manipulated during scaling and then used to perform extended precision
     * computation of the sum {@code x^2 + y^2} where the high part of the number has 20-bit
     * precision. Manipulation of direct bits has no equivalent in Java
     * other than use of {@link Double#doubleToLongBits(double)} and
     * {@link Double#longBitsToDouble(long)}. To avoid conversion to and from long and double
     * representations this implementation only scales the double representation. The high
     * and low parts of a double for the extended precision computation are extracted
     * using the method of Dekker (1971) to create two 26-bit numbers. This works for sub-normal
     * numbers and reduces the maximum error in comparison to fdlibm hypot which does not
     * use a split number algorithm for sub-normal numbers.
     *
     * @param x Value x
     * @param y Value y
     * @return sqrt(x^2 + y^2)
     * @see Math#hypot(double, double)
     * @see <a href="https://www.netlib.org/fdlibm/e_hypot.c">fdlibm e_hypot.c</a>
     * @see <a href="https://bugs.java.com/bugdatabase/view_bug.do?bug_id=7130085">JDK-7130085 : Port fdlibm hypot to Java</a>
     */
    public static double hypot(double x, double y) {
        // Differences to the fdlibm reference:
        //
        // 1. fdlibm orders the two parts using the magnitude of the upper 32-bits.
        // This incorrectly orders numbers which differ only in the lower 32-bits.
        // This invalidates hypot(x, y) = hypot(y, x) for small sub-normal numbers and a minority
        // of cases of normal numbers. This implementation forces the |x| >= |y| order
        // using the entire 63-bits of the unsigned doubles to ensure the function
        // is commutative.
        //
        // 2. fdlibm computed scaling by directly writing changes to the exponent bits
        // and maintained the high part (ha) during scaling for use in the high
        // precision sum x^2 + y^2. Since exponent scaling cannot be applied to sub-normals
        // the original version dropped the split number representation for sub-normals
        // and can produce maximum errors above 1 ULP for sub-normal numbers.
        // This version uses Dekker's method to split the number. This can be applied to
        // sub-normals and allows dropping the condition to check for sub-normal numbers
        // since all small numbers are handled with a single scaling factor.
        // The effect is increased precision for the majority of sub-normal cases where
        // the implementations compute a different result.
        //
        // 3. An alteration is done here to add an 'else if' instead of a second
        // 'if' statement. Thus you cannot scale down and up at the same time.
        //
        // 4. There is no use of the absolute double value. The magnitude comparison is
        // performed using the long bit representation. The computation x^2+y^2 is
        // insensitive to the sign bit. Thus use of Math.abs(double) is only in edge-case
        // branches.
        //
        // 5. The exponent different to ignore the smaller component has changed from 60 to 54.
        //
        // Original comments from fdlibm are in c style: /* */
        // Extra comments added for reference.
        //
        // Note that the high 32-bits are compared to constants.
        // The lowest 20-bits are the upper bits of the 52-bit mantissa.
        // The next 11-bits are the biased exponent. The sign bit has been cleared.
        // Scaling factors are powers of two for exact scaling.
        // For clarity the values have been refactored to named constants.

        // The mask is used to remove the sign bit.
        final long xbits = Double.doubleToRawLongBits(x) & UNSIGN_MASK;
        final long ybits = Double.doubleToRawLongBits(y) & UNSIGN_MASK;

        // Order by magnitude: |a| >= |b|
        double a;
        double b;
        /* High word of x & y */
        int ha;
        int hb;
        if (ybits > xbits) {
            a = y;
            b = x;
            ha = (int) (ybits >>> 32);
            hb = (int) (xbits >>> 32);
        } else {
            a = x;
            b = y;
            ha = (int) (xbits >>> 32);
            hb = (int) (ybits >>> 32);
        }

        // Check if the smaller part is significant.
        // a^2 is computed in extended precision for an effective mantissa of 106-bits.
        // An exponent difference of 54 is where b^2 will not overlap a^2.
        if ((ha - hb) > EXP_54) {
            /* a/b > 2**54 */
            // or a is Inf or NaN.
            // No addition of a + b for sNaN.
            return Math.abs(a);
        }

        double rescale = 1.0;
        if (ha > EXP_500) {
            /* a > 2^500 */
            if (ha >= EXP_1024) {
                /* Inf or NaN */
                // Check b is infinite for the IEEE754 result.
                // No addition of a + b for sNaN.
                return Math.abs(b) == Double.POSITIVE_INFINITY ?
                    Double.POSITIVE_INFINITY :
                    Math.abs(a);
            }
            /* scale a and b by 2^-600 */
            // Before scaling: a in [2^500, 2^1023].
            // After scaling: a in [2^-100, 2^423].
            // After scaling: b in [2^-154, 2^423].
            a *= TWO_POW_NEG_600;
            b *= TWO_POW_NEG_600;
            rescale = TWO_POW_600;
        } else if (hb < EXP_NEG_500) {
            // No special handling of sub-normals.
            // These do not matter when we do not manipulate the exponent bits
            // for scaling the split representation.

            // Intentional comparison with zero.
            if (b == 0) {
                return Math.abs(a);
            }

            /* scale a and b by 2^600 */
            // Effective min exponent of a sub-normal = -1022 - 52 = -1074.
            // Before scaling: b in [2^-1074, 2^-501].
            // After scaling: b in [2^-474, 2^99].
            // After scaling: a in [2^-474, 2^153].
            a *= TWO_POW_600;
            b *= TWO_POW_600;
            rescale = TWO_POW_NEG_600;
        }

        // High precision x^2 + y^2
        return Math.sqrt(x2y2(a, b)) * rescale;
    }

    private static double x2y2(double x, double y) {
        // Note:
        // This method is different from the high-accuracy summation used in fdlibm for hypot.
        // The summation could be any valid computation of x^2+y^2. However since this follows
        // the re-scaling logic in hypot(x, y) the use of high precision has relatively
        // less performance overhead than if used without scaling.
        // The Dekker algorithm is branchless for better performance
        // than the fdlibm method with a maximum ULP error of approximately 0.86.
        //
        // See NUMBERS-143 for analysis.

        // Do a Dekker summation of double length products x*x and y*y
        // (10 multiply and 20 additions).
        final double xx = x * x;
        final double yy = y * y;
        // Compute the round-off from the products.
        // With FMA hardware support in JDK 9+ this can be replaced with the much faster:
        // xxLow = Math.fma(x, x, -xx)
        // yyLow = Math.fma(y, y, -yy)
        // Dekker mul12
        final double xHigh = splitHigh(x);
        final double xLow = x - xHigh;
        final double xxLow = squareLow(xLow, xHigh, xx);
        // Dekker mul12
        final double yHigh = splitHigh(y);
        final double yLow = y - yHigh;
        final double yyLow = squareLow(yLow, yHigh, yy);
        // Dekker add2
        final double r = xx + yy;
        // Note: The order is important. Assume xx > yy and drop Dekker's conditional
        // check for which is the greater magnitude.
        // s = xx - r + yy + yyLow + xxLow
        // z = r + s
        // zz = r - z + s
        // Here we compute z inline and ignore computing the round-off zz.
        // Note: The round-off could be used with Dekker's sqrt2 method.
        // That adds 7 multiply, 1 division and 19 additions doubling the cost
        // and reducing error to < 0.5 ulp for the final sqrt.
        return xx - r + yy + yyLow + xxLow + r;
    }
    public static <R> R add(double r1, double i1, double r2, double i2, ComplexResult<R> result) {
        return result.apply(r1 + r2, i1 + i2);
    }

    /**
     * Returns a {@code Complex} whose value is:
     * <pre>
     *  (a + i b)(c + i d) = (ac - bd) + i (ad + bc)</pre>
     *
     * <p>Recalculates to recover infinities as specified in C99 standard G.5.1.
     *
     * @param r1 Real component of first number.
     * @param i1 Imaginary component of first number.
     * @param r2 Real component of second number.
     * @param i2 Imaginary component of second number.
     * @param <R> Generic
     * @param result Complex result
     * @return (a + b i)(c + d i).
     */
    public static <R> R multiply(double r1, double i1, double r2, double i2, ComplexResult<R> result) {
        double a = r1;
        double b = i1;
        double c = r2;
        double d = i2;
        final double ac = a * c;
        final double bd = b * d;
        final double ad = a * d;
        final double bc = b * c;
        double x = ac - bd;
        double y = ad + bc;

        // --------------
        // NaN can occur if:
        // - any of (a,b,c,d) are NaN (for NaN or Infinite complex numbers)
        // - a multiplication of infinity by zero (ac,bd,ad,bc).
        // - a subtraction of infinity from infinity (e.g. ac - bd)
        //   Note that (ac,bd,ad,bc) can be infinite due to overflow.
        //
        // Detect a NaN result and perform correction.
        //
        // Modification from the listing in ISO C99 G.5.1 (6)
        // Do not correct infinity multiplied by zero. This is left as NaN.
        // --------------

        if (Double.isNaN(x) && Double.isNaN(y)) {
            // Recover infinities that computed as NaN+iNaN ...
            boolean recalc = false;
            if ((Double.isInfinite(a) || Double.isInfinite(b)) &&
                isNotZero(c, d)) {
                // This complex is infinite.
                // "Box" the infinity and change NaNs in the other factor to 0.
                a = boxInfinity(a);
                b = boxInfinity(b);
                c = changeNaNtoZero(c);
                d = changeNaNtoZero(d);
                recalc = true;
            }
            if ((Double.isInfinite(c) || Double.isInfinite(d)) &&
                isNotZero(a, b)) {
                // The other complex is infinite.
                // "Box" the infinity and change NaNs in the other factor to 0.
                c = boxInfinity(c);
                d = boxInfinity(d);
                a = changeNaNtoZero(a);
                b = changeNaNtoZero(b);
                recalc = true;
            }
            if (!recalc && (Double.isInfinite(ac) || Double.isInfinite(bd) ||
                Double.isInfinite(ad) || Double.isInfinite(bc))) {
                // The result overflowed to infinity.
                // Recover infinities from overflow by changing NaNs to 0 ...
                a = changeNaNtoZero(a);
                b = changeNaNtoZero(b);
                c = changeNaNtoZero(c);
                d = changeNaNtoZero(d);
                recalc = true;
            }
            if (recalc) {
                x = Double.POSITIVE_INFINITY * (a * c - b * d);
                y = Double.POSITIVE_INFINITY * (a * d + b * c);
            }
        }
        return result.apply(x, y);
    }
    /**
     * Box values for the real or imaginary component of an infinite complex number.
     * Any infinite value will be returned as one. Non-infinite values will be returned as zero.
     * The sign is maintained.
     *
     * <pre>
     *  inf  =  1
     * -inf  = -1
     *  x    =  0
     * -x    = -0
     * </pre>
     *
     * @param component the component
     * @return The boxed value
     */
    private static double boxInfinity(double component) {
        return Math.copySign(Double.isInfinite(component) ? 1.0 : 0.0, component);
    }

    /**
     * Checks if the complex number is not zero.
     *
     * @param real the real component
     * @param imaginary the imaginary component
     * @return true if the complex is not zero
     */
    private static boolean isNotZero(double real, double imaginary) {
        // The use of equals is deliberate.
        // This method must distinguish NaN from zero thus ruling out:
        // (real != 0.0 || imaginary != 0.0)
        return !(real == 0.0 && imaginary == 0.0);
    }

    /**
     * Change NaN to zero preserving the sign; otherwise return the value.
     *
     * @param value the value
     * @return The new value
     */
    private static double changeNaNtoZero(double value) {
        return Double.isNaN(value) ? Math.copySign(0.0, value) : value;
    }


    /**
     * Returns a {@code Complex} whose value is:
     * <pre>
     * <code>
     *   a + i b     (ac + bd) + i (bc - ad)
     *   -------  =  -----------------------
     *   c + i d            c<sup>2</sup> + d<sup>2</sup>
     * </code>
     * </pre>
     *
     * <p>Recalculates to recover infinities as specified in C99
     * standard G.5.1. Method is fully in accordance with
     * C++11 standards for complex numbers.</p>
     *
     * <p>Note: In the event of divide by zero this method produces the same result
     * as dividing by a real-only zero using {@link #(double)}.
     *
     * @param re1 Real component of first number.
     * @param im1 Imaginary component of first number.
     * @param re2 Real component of second number.
     * @param im2 Imaginary component of second number.
     * @param <R> Generic
     * @param result Complex result
     * @return (a + i b) / (c + i d).
     * @see <a href="http://mathworld.wolfram.com/ComplexDivision.html">Complex Division</a>
     * @see #(double)
     */
    public static <R> R divide(double re1, double im1, double re2, double im2, ComplexResult<R> result) {
        double a = re1;
        double b = im1;
        double c = re2;
        double d = im2;
        int ilogbw = 0;
        // Get the exponent to scale the divisor parts to the range [1, 2).
        final int exponent = getScale(c, d);
        if (exponent <= Double.MAX_EXPONENT) {
            ilogbw = exponent;
            c = Math.scalb(c, -ilogbw);
            d = Math.scalb(d, -ilogbw);
        }
        final double denom = c * c + d * d;

        // Note: Modification from the listing in ISO C99 G.5.1 (8):
        // Avoid overflow if a or b are very big.
        // Since (c, d) in the range [1, 2) the sum (ac + bd) could overflow
        // when (a, b) are both above (Double.MAX_VALUE / 4). The same applies to
        // (bc - ad) with large negative values.
        // Use the maximum exponent as an approximation to the magnitude.
        if (getMaxExponent(a, b) > Double.MAX_EXPONENT - 2) {
            ilogbw -= 2;
            a /= 4;
            b /= 4;
        }

        double x = Math.scalb((a * c + b * d) / denom, -ilogbw);
        double y = Math.scalb((b * c - a * d) / denom, -ilogbw);
        // Recover infinities and zeros that computed as NaN+iNaN
        // the only cases are nonzero/zero, infinite/finite, and finite/infinite, ...
        if (Double.isNaN(x) && Double.isNaN(y)) {
            if ((denom == 0.0) &&
                (!Double.isNaN(a) || !Double.isNaN(b))) {
                // nonzero/zero
                // This case produces the same result as divide by a real-only zero
                // using Complex.divide(+/-0.0)
                x = Math.copySign(Double.POSITIVE_INFINITY, c) * a;
                y = Math.copySign(Double.POSITIVE_INFINITY, c) * b;
            } else if ((Double.isInfinite(a) || Double.isInfinite(b)) &&
                Double.isFinite(c) && Double.isFinite(d)) {
                // infinite/finite
                a = boxInfinity(a);
                b = boxInfinity(b);
                x = Double.POSITIVE_INFINITY * (a * c + b * d);
                y = Double.POSITIVE_INFINITY * (b * c - a * d);
            } else if ((Double.isInfinite(c) || Double.isInfinite(d)) &&
                Double.isFinite(a) && Double.isFinite(b)) {
                // finite/infinite
                c = boxInfinity(c);
                d = boxInfinity(d);
                x = 0.0 * (a * c + b * d);
                y = 0.0 * (b * c - a * d);
            }
        }
        return result.apply(x, y);
    }
    /**
     * Returns the
     * <a href="http://mathworld.wolfram.com/ExponentialFunction.html">
     * exponential function</a> of this complex number.
     *
     * <p>\[ \exp(z) = e^z \]
     *
     * <p>The exponential function of \( z \) is an entire function in the complex plane.
     * Special cases:
     *
     * <ul>
     * <li>{@code z.conj().exp() == z.exp().conj()}.
     * <li>If {@code z} is ±0 + i0, returns 1 + i0.
     * <li>If {@code z} is x + i∞ for finite x, returns NaN + iNaN ("invalid" floating-point operation).
     * <li>If {@code z} is x + iNaN for finite x, returns NaN + iNaN ("invalid" floating-point operation).
     * <li>If {@code z} is +∞ + i0, returns +∞ + i0.
     * <li>If {@code z} is −∞ + iy for finite y, returns +0 cis(y).
     * <li>If {@code z} is +∞ + iy for finite nonzero y, returns +∞ cis(y).
     * <li>If {@code z} is −∞ + i∞, returns ±0 ± i0 (where the signs of the real and imaginary parts of the result are unspecified).
     * <li>If {@code z} is +∞ + i∞, returns ±∞ + iNaN (where the sign of the real part of the result is unspecified; "invalid" floating-point operation).
     * <li>If {@code z} is −∞ + iNaN, returns ±0 ± i0 (where the signs of the real and imaginary parts of the result are unspecified).
     * <li>If {@code z} is +∞ + iNaN, returns ±∞ + iNaN (where the sign of the real part of the result is unspecified).
     * <li>If {@code z} is NaN + i0, returns NaN + i0.
     * <li>If {@code z} is NaN + iy for all nonzero numbers y, returns NaN + iNaN ("invalid" floating-point operation).
     * <li>If {@code z} is NaN + iNaN, returns NaN + iNaN.
     * </ul>
     *
     * <p>Implements the formula:
     *
     * <p>\[ \exp(x + iy) = e^x (\cos(y) + i \sin(y)) \]
     * @param <R> Generic
     * @param result Complex result
     * @param r Real part
     * @param i Imaginary part
     * @return The exponential of this complex number.
     * @see <a href="http://functions.wolfram.com/ElementaryFunctions/Exp/">Exp</a>
     */

    public static <R> R exp(double r, double i, ComplexResult<R> result) {
        if (Double.isInfinite(r)) {
            // Set the scale factor applied to cis(y)
            double zeroOrInf;
            if (r < 0) {
                if (!Double.isFinite(i)) {
                    // (−∞ + i∞) or (−∞ + iNaN) returns (±0 ± i0) (where the signs of the
                    // real and imaginary parts of the result are unspecified).
                    // Here we preserve the conjugate equality.
                    return result.apply(0, Math.copySign(0, i));
                }
                // (−∞ + iy) returns +0 cis(y), for finite y
                zeroOrInf = 0;
            } else {
                // (+∞ + i0) returns +∞ + i0.
                if (i == 0) {
                    return result.apply(r, i);
                }
                // (+∞ + i∞) or (+∞ + iNaN) returns (±∞ + iNaN) and raises the invalid
                // floating-point exception (where the sign of the real part of the
                // result is unspecified).
                if (!Double.isFinite(i)) {
                    return result.apply(r, Double.NaN);
                }
                // (+∞ + iy) returns (+∞ cis(y)), for finite nonzero y.
                zeroOrInf = r;
            }
            return result.apply(zeroOrInf * Math.cos(i),
                zeroOrInf * Math.sin(i));

        } else if (Double.isNaN(r)) {
            // (NaN + i0) returns (NaN + i0)
            // (NaN + iy) returns (NaN + iNaN) and optionally raises the invalid floating-point exception
            // (NaN + iNaN) returns (NaN + iNaN)
            if (i == 0) {
                return result.apply(r, i);
            } else {
                return result.apply(Double.NaN, Double.NaN);
            }
        } else if (!Double.isFinite(i)) {
            // (x + i∞) or (x + iNaN) returns (NaN + iNaN) and raises the invalid
            // floating-point exception, for finite x.
            return result.apply(Double.NaN, Double.NaN);
        }
        // real and imaginary are finite.
        // Compute e^a * (cos(b) + i sin(b)).

        // Special case:
        // (±0 + i0) returns (1 + i0)
        final double exp = Math.exp(r);
        if (i == 0) {
            return result.apply(exp, i);
        }
        return result.apply(exp * Math.cos(i),
            exp * Math.sin(i));

    }

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

    /**
     * Returns the inverse sine of the complex number.
     *
     * <p>This function exists to allow implementation of the identity
     * {@code asinh(z) = -i asin(iz)}.
     *
     * <p>Adapted from {@code <boost/math/complex/asin.hpp>}. This method only (and not
     * invoked methods within) is distributed under the Boost Software License V1.0.
     * The original notice is shown below and the licence is shown in full in LICENSE:
     * <pre>
     * (C) Copyright John Maddock 2005.
     * Distributed under the Boost Software License, Version 1.0. (See accompanying
     * file LICENSE or copy at https://www.boost.org/LICENSE_1_0.txt)
     * </pre>
     *
     * @param real Real part.
     * @param imaginary Imaginary part.
     * @param result Complex Result.
     * @param <R> Generic
     * @return The inverse sine of this complex number.
     */
    public static <R> R asin(final double real, final double imaginary, ComplexResult<R> result) {
        // Compute with positive values and determine sign at the end
        final double x = Math.abs(real);
        final double y = Math.abs(imaginary);
        // The result (without sign correction)
        double re;
        double im;

        // Handle C99 special cases
        if (Double.isNaN(x)) {
            if (isPosInfinite(y)) {
                re = x;
                im = y;
            } else {
                // No-use of the input constructor
                return result.apply(Double.NaN, Double.NaN);
            }
        } else if (Double.isNaN(y)) {
            if (x == 0) {
                re = 0;
                im = y;
            } else if (isPosInfinite(x)) {
                re = y;
                im = x;
            } else {
                // No-use of the input constructor
                return result.apply(Double.NaN, Double.NaN);
            }
        } else if (isPosInfinite(x)) {
            re = isPosInfinite(y) ? PI_OVER_4 : PI_OVER_2;
            im = x;
        } else if (isPosInfinite(y)) {
            re = 0;
            im = y;
        } else {
            // Special case for real numbers:
            if (y == 0 && x <= 1) {
                return result.apply(Math.asin(real), imaginary);
            }

            final double xp1 = x + 1;
            final double xm1 = x - 1;

            if (inRegion(x, y, SAFE_MIN, SAFE_MAX)) {
                final double yy = y * y;
                final double r = Math.sqrt(xp1 * xp1 + yy);
                final double s = Math.sqrt(xm1 * xm1 + yy);
                final double a = 0.5 * (r + s);
                final double b = x / a;

                if (b <= B_CROSSOVER) {
                    re = Math.asin(b);
                } else {
                    final double apx = a + x;
                    if (x <= 1) {
                        re = Math.atan(x / Math.sqrt(0.5 * apx * (yy / (real + xp1) + (s - xm1))));
                    } else {
                        re = Math.atan(x / (y * Math.sqrt(0.5 * (apx / (real + xp1) + apx / (s + xm1)))));
                    }
                }

                if (a <= A_CROSSOVER) {
                    double am1;
                    if (x < 1) {
                        am1 = 0.5 * (yy / (real + xp1) + yy / (s - xm1));
                    } else {
                        am1 = 0.5 * (yy / (real + xp1) + (s + xm1));
                    }
                    im = Math.log1p(am1 + Math.sqrt(am1 * (a + 1)));
                } else {
                    im = Math.log(a + Math.sqrt(a * a - 1));
                }
            } else {
                // Hull et al: Exception handling code from figure 4
                if (y <= (EPSILON * Math.abs(xm1))) {
                    if (x < 1) {
                        re = Math.asin(x);
                        im = y / Math.sqrt(xp1 * (1 - x));
                    } else {
                        re = PI_OVER_2;
                        if ((Double.MAX_VALUE / xp1) > xm1) {
                            // xp1 * xm1 won't overflow:
                            im = Math.log1p(xm1 + Math.sqrt(xp1 * xm1));
                        } else {
                            im = LN_2 + Math.log(x);
                        }
                    }
                } else if (y <= SAFE_MIN) {
                    // Hull et al: Assume x == 1.
                    // True if:
                    // E^2 > 8*sqrt(u)
                    //
                    // E = Machine epsilon: (1 + epsilon) = 1
                    // u = Double.MIN_NORMAL
                    re = PI_OVER_2 - Math.sqrt(y);
                    im = Math.sqrt(y);
                } else if (EPSILON * y - 1 >= x) {
                    // Possible underflow:
                    re = x / y;
                    im = LN_2 + Math.log(y);
                } else if (x > 1) {
                    re = Math.atan(x / y);
                    final double xoy = x / y;
                    im = LN_2 + Math.log(y) + 0.5 * Math.log1p(xoy * xoy);
                } else {
                    final double a = Math.sqrt(1 + y * y);
                    // Possible underflow:
                    re = x / a;
                    im = 0.5 * Math.log1p(2 * y * (y + a));
                }
            }
        }
        return result.apply(changeSign(re, real),
            changeSign(im, imaginary));
    }

    /**
     * Returns the inverse cosine of the complex number.
     *
     * <p>This function exists to allow implementation of the identity
     * {@code acosh(z) = +-i acos(z)}.
     *
     * <p>Adapted from {@code <boost/math/complex/acos.hpp>}. This method only (and not
     * invoked methods within) is distributed under the Boost Software License V1.0.
     * The original notice is shown below and the licence is shown in full in LICENSE:
     * <pre>
     * (C) Copyright John Maddock 2005.
     * Distributed under the Boost Software License, Version 1.0. (See accompanying
     * file LICENSE or copy at https://www.boost.org/LICENSE_1_0.txt)
     * </pre>
     *
     * @param real Real part.
     * @param imaginary Imaginary part.
     * @param result Constructor.
     * @return The inverse cosine of the complex number.
     */
    public static <R> R acos(final double real, final double imaginary,
                                final ComplexResult<R> result) {
        // Compute with positive values and determine sign at the end
        final double x = Math.abs(real);
        final double y = Math.abs(imaginary);
        // The result (without sign correction)
        double re;
        double im;

        // Handle C99 special cases
        if (isPosInfinite(x)) {
            if (isPosInfinite(y)) {
                re = PI_OVER_4;
                im = y;
            } else if (Double.isNaN(y)) {
                // sign of the imaginary part of the result is unspecified
                return result.apply(imaginary, real);
            } else {
                re = 0;
                im = Double.POSITIVE_INFINITY;
            }
        } else if (Double.isNaN(x)) {
            if (isPosInfinite(y)) {
                return result.apply(x, -imaginary);
            }
            // No-use of the input constructor
            return (R) NAN;
        } else if (isPosInfinite(y)) {
            re = PI_OVER_2;
            im = y;
        } else if (Double.isNaN(y)) {
            return result.apply(x == 0 ? PI_OVER_2 : y, y);
        } else {
            // Special case for real numbers:
            if (y == 0 && x <= 1) {
                return result.apply(x == 0 ? PI_OVER_2 : Math.acos(real), -imaginary);
            }

            final double xp1 = x + 1;
            final double xm1 = x - 1;

            if (inRegion(x, y, SAFE_MIN, SAFE_MAX)) {
                final double yy = y * y;
                final double r = Math.sqrt(xp1 * xp1 + yy);
                final double s = Math.sqrt(xm1 * xm1 + yy);
                final double a = 0.5 * (r + s);
                final double b = x / a;

                if (b <= B_CROSSOVER) {
                    re = Math.acos(b);
                } else {
                    final double apx = a + x;
                    if (x <= 1) {
                        re = Math.atan(Math.sqrt(0.5 * apx * (yy / (r + xp1) + (s - xm1))) / x);
                    } else {
                        re = Math.atan((y * Math.sqrt(0.5 * (apx / (r + xp1) + apx / (s + xm1)))) / x);
                    }
                }

                if (a <= A_CROSSOVER) {
                    double am1;
                    if (x < 1) {
                        am1 = 0.5 * (yy / (r + xp1) + yy / (s - xm1));
                    } else {
                        am1 = 0.5 * (yy / (r + xp1) + (s + xm1));
                    }
                    im = Math.log1p(am1 + Math.sqrt(am1 * (a + 1)));
                } else {
                    im = Math.log(a + Math.sqrt(a * a - 1));
                }
            } else {
                // Hull et al: Exception handling code from figure 6
                if (y <= (EPSILON * Math.abs(xm1))) {
                    if (x < 1) {
                        re = Math.acos(x);
                        im = y / Math.sqrt(xp1 * (1 - x));
                    } else {
                        // This deviates from Hull et al's paper as per
                        // https://svn.boost.org/trac/boost/ticket/7290
                        if ((Double.MAX_VALUE / xp1) > xm1) {
                            // xp1 * xm1 won't overflow:
                            re = y / Math.sqrt(xm1 * xp1);
                            im = Math.log1p(xm1 + Math.sqrt(xp1 * xm1));
                        } else {
                            re = y / x;
                            im = LN_2 + Math.log(x);
                        }
                    }
                } else if (y <= SAFE_MIN) {
                    // Hull et al: Assume x == 1.
                    // True if:
                    // E^2 > 8*sqrt(u)
                    //
                    // E = Machine epsilon: (1 + epsilon) = 1
                    // u = Double.MIN_NORMAL
                    re = Math.sqrt(y);
                    im = Math.sqrt(y);
                } else if (EPSILON * y - 1 >= x) {
                    re = PI_OVER_2;
                    im = LN_2 + Math.log(y);
                } else if (x > 1) {
                    re = Math.atan(y / x);
                    final double xoy = x / y;
                    im = LN_2 + Math.log(y) + 0.5 * Math.log1p(xoy * xoy);
                } else {
                    re = PI_OVER_2;
                    final double a = Math.sqrt(1 + y * y);
                    im = 0.5 * Math.log1p(2 * y * (y + a));
                }
            }
        }

        return result.apply(negative(real) ? Math.PI - re : re,
            negative(imaginary) ? im : -im);
    }
    public static <R> R acosh(double r, double i, ComplexResult<R> result) {
        // Define in terms of acos
        // acosh(z) = +-i acos(z)
        // Note the special case:
        // acos(+-0 + iNaN) = π/2 + iNaN
        // acosh(0 + iNaN) = NaN + iπ/2
        // will not appropriately multiply by I to maintain positive imaginary if
        // acos() imaginary computes as NaN. So do this explicitly.
        if (Double.isNaN(i) && r == 0) {
            return result.apply(Double.NaN, PI_OVER_2);
        }
        return acos(r, i, (re, im) ->
            // Set the sign appropriately for real >= 0
            (negative(im)) ?
                // Multiply by I
                result.apply(-im, re) :
                // Multiply by -I
                result.apply(im, -re)
        );
    }
    public static <R> R conj(double r, double i, ComplexResult<R> result) {
        return result.apply(r, -i);
    }

    public static <R> R add(double r, double i, double a, ComplexResult<R> result) {
        return result.apply(r + a, i);
    }
    public static <R> R addImaginary(double r, double i, double a, ComplexResult<R> result) {
        return result.apply(r, i + a);
    }
    public static <R> R subtract(double r, double i, double s, ComplexResult<R> result) {
        return result.apply(r - s, i);
    }
    public static <R> R subtractImaginary(double r, double i, double s, ComplexResult<R> result) {
        return result.apply(r, i - s);
    }
    public static <R> R subtractFrom(double r, double i, double m, ComplexResult<R> result) {
        return result.apply(m - r, -i);
    }
    public static <R> R subtractFromImaginary(double r, double i, double m, ComplexResult<R> result) {
        return result.apply(-r, m - i);
    }
    public static <R> R multiply(double r, double i, double f, ComplexResult<R> result)
    {
        return result.apply(r * f, i * f);
    }
    public static <R> R multiplyImaginary(double r, double i, double f, ComplexResult<R> result) {
        return result.apply(-i * f, r * f);
    }
    public static <R> R divide(double r, double i, double d, ComplexResult<R> result) {
        return result.apply(r / d, i / d);
    }
    public static <R> R divideImaginary(double r, double i, double d, ComplexResult<R> result) {
        return result.apply(i / d, -r / d);
    }
    /**
     * Returns the logarithm of this complex number using the provided function.
     * Implements the formula:
     *
     * <pre>
     *   log(x + i y) = log(|x + i y|) + i arg(x + i y)</pre>
     *
     * <p>Warning: The argument {@code logOf2} must be equal to {@code log(2)} using the
     * provided log function otherwise scaling using powers of 2 in the case of overflow
     * will be incorrect. This is provided as an internal optimisation.
     *
     * @param log Log function.
     * @param logOfeOver2 The log function applied to e, then divided by 2.
     * @param logOf2 The log function applied to 2.
     * @param result Constructor for the returned complex.
     * @return The logarithm of this complex number.
     * @see #abs(double, double)
     * @see #arg(double, double)
     */
    private static <R> R log(double real, double imaginary, DoubleUnaryOperator log, double logOfeOver2, double logOf2, ComplexResult<R> result) {
        // Handle NaN
        if (Double.isNaN(real) || Double.isNaN(imaginary)) {
            // Return NaN unless infinite
            if (isInfinite(real, imaginary)) {
                return result.apply(Double.POSITIVE_INFINITY, Double.NaN);
            }
            // No-use of the input constructor
            return (R) NAN;
        }

        // Returns the real part:
        // log(sqrt(x^2 + y^2))
        // log(x^2 + y^2) / 2

        // Compute with positive values
        double x = Math.abs(real);
        double y = Math.abs(imaginary);

        // Find the larger magnitude.
        if (x < y) {
            final double tmp = x;
            x = y;
            y = tmp;
        }

        if (x == 0) {
            // Handle zero: raises the ‘‘divide-by-zero’’ floating-point exception.
            return result.apply(Double.NEGATIVE_INFINITY,
                negative(real) ? Math.copySign(Math.PI, imaginary) : imaginary);
        }

        double re;

        // This alters the implementation of Hull et al (1994) which used a standard
        // precision representation of |z|: sqrt(x*x + y*y).
        // This formula should use the same definition of the magnitude returned
        // by Complex.abs() which is a high precision computation with scaling.
        // The checks for overflow thus only require ensuring the output of |z|
        // will not overflow or underflow.

        if (x > HALF && x < ROOT2) {
            // x^2+y^2 close to 1. Use log1p(x^2+y^2 - 1) / 2.
            re = Math.log1p(x2y2m1(x, y)) * logOfeOver2;
        } else {
            // Check for over/underflow in |z|
            // When scaling:
            // log(a / b) = log(a) - log(b)
            // So initialize the result with the log of the scale factor.
            re = 0;
            if (x > Double.MAX_VALUE / 2) {
                // Potential overflow.
                if (isPosInfinite(x)) {
                    // Handle infinity
                    return result.apply(x, arg(real, imaginary));
                }
                // Scale down.
                x /= 2;
                y /= 2;
                // log(2)
                re = logOf2;
            } else if (y < Double.MIN_NORMAL) {
                // Potential underflow.
                if (y == 0) {
                    // Handle real only number
                    return result.apply(log.applyAsDouble(x), arg(real, imaginary));
                }
                // Scale up sub-normal numbers to make them normal by scaling by 2^54,
                // i.e. more than the mantissa digits.
                x *= 0x1.0p54;
                y *= 0x1.0p54;
                // log(2^-54) = -54 * log(2)
                re = -54 * logOf2;
            }
            re += log.applyAsDouble(abs(x, y));
        }

        // All ISO C99 edge cases for the imaginary are satisfied by the Math library.
        return result.apply(re, arg(real, imaginary));
    }
    public static <R> R sqrt(double real, double imaginary, ComplexResult<R> result) {
        // Handle NaN
        if (Double.isNaN(real) || Double.isNaN(imaginary)) {
            // Check for infinite
            if (Double.isInfinite(imaginary)) {
                return result.apply(Double.POSITIVE_INFINITY, imaginary);
            }
            if (Double.isInfinite(real)) {
                if (real == Double.NEGATIVE_INFINITY) {
                    return result.apply(Double.NaN, Math.copySign(Double.POSITIVE_INFINITY, imaginary));
                }
                return result.apply(Double.POSITIVE_INFINITY, Double.NaN);
            }
            return result.apply(Double.NaN, Double.NaN);
        }

        // Compute with positive values and determine sign at the end
        final double x = Math.abs(real);
        final double y = Math.abs(imaginary);

        // Compute
        double t;

        // This alters the implementation of Hull et al (1994) which used a standard
        // precision representation of |z|: sqrt(x*x + y*y).
        // This formula should use the same definition of the magnitude returned
        // by Complex.abs() which is a high precision computation with scaling.
        // Worry about overflow if 2 * (|z| + |x|) will overflow.
        // Worry about underflow if |z| or |x| are sub-normal components.

        if (inRegion(x, y, Double.MIN_NORMAL, SQRT_SAFE_UPPER)) {
            // No over/underflow
            t = Math.sqrt(2 * (abs(x, y) + x));
        } else {
            // Potential over/underflow. First check infinites and real/imaginary only.

            // Check for infinite
            if (isPosInfinite(y)) {
                return result.apply(Double.POSITIVE_INFINITY, imaginary);
            } else if (isPosInfinite(x)) {
                if (real == Double.NEGATIVE_INFINITY) {
                    return result.apply(0, Math.copySign(Double.POSITIVE_INFINITY, imaginary));
                }
                return result.apply(Double.POSITIVE_INFINITY, Math.copySign(0, imaginary));
            } else if (y == 0) {
                // Real only
                final double sqrtAbs = Math.sqrt(x);
                if (real < 0) {
                    return result.apply(0, Math.copySign(sqrtAbs, imaginary));
                }
                return result.apply(sqrtAbs, imaginary);
            } else if (x == 0) {
                // Imaginary only. This sets the two components to the same magnitude.
                // Note: In polar coordinates this does not happen:
                // real = sqrt(abs()) * Math.cos(arg() / 2)
                // imag = sqrt(abs()) * Math.sin(arg() / 2)
                // arg() / 2 = pi/4 and cos and sin should both return sqrt(2)/2 but
                // are different by 1 ULP.
                final double sqrtAbs = Math.sqrt(y) * ONE_OVER_ROOT2;
                return result.apply(sqrtAbs, Math.copySign(sqrtAbs, imaginary));
            } else {
                // Over/underflow.
                // Full scaling is not required as this is done in the hypotenuse function.
                // Keep the number as big as possible for maximum precision in the second sqrt.
                // Note if we scale by an even power of 2, we can re-scale by sqrt of the number.
                // a = sqrt(b)
                // a = sqrt(b/4) * sqrt(4)

                double rescale;
                double sx;
                double sy;
                if (Math.max(x, y) > SQRT_SAFE_UPPER) {
                    // Overflow. Scale down by 16 and rescale by sqrt(16).
                    sx = x / 16;
                    sy = y / 16;
                    rescale = 4;
                } else {
                    // Sub-normal numbers. Make them normal by scaling by 2^54,
                    // i.e. more than the mantissa digits, and rescale by sqrt(2^54) = 2^27.
                    sx = x * 0x1.0p54;
                    sy = y * 0x1.0p54;
                    rescale = 0x1.0p-27;
                }
                t = rescale * Math.sqrt(2 * (abs(sx, sy) + sx));
            }
        }

        if (real >= 0) {
            return result.apply(t / 2, imaginary / t);
        }
        return result.apply(y / t, Math.copySign(t / 2, imaginary));
    }

    /*
    public static <R> R pow(double r1, double i1, double r2, double i2, ComplexResult<R> result) {
        if (r1 == 0 &&
            i1 == 0) {
            // This value is zero. Test the other.
            if (r2 > 0 &&
                i2 == 0) {
                // 0 raised to positive number is 0
                return result.apply(0, 0);
            }
            // 0 raised to anything else is NaN
            return result.apply(Double.NaN, Double.NaN);
        }
        return log().multiply().exp();
    }
    */



    public static <R> R subtract(double re1, double im1, double re2, double im2, ComplexResult<R> result) {
        return result.apply(re1 - re2,
            im1 - im2);
    }

    public static <R> R sinh(double real, double imaginary, ComplexResult<R> result) {
        if (Double.isInfinite(real) && !Double.isFinite(imaginary)) {
            return result.apply(real, Double.NaN);
        }
        if (real == 0) {
            // Imaginary-only sinh(iy) = i sin(y).
            if (Double.isFinite(imaginary)) {
                // Maintain periodic property with respect to the imaginary component.
                // sinh(+/-0.0) * cos(+/-x) = +/-0 * cos(x)
                return result.apply(changeSign(real, Math.cos(imaginary)),
                    Math.sin(imaginary));
            }
            // If imaginary is inf/NaN the sign of the real part is unspecified.
            // Returning the same real value maintains the conjugate equality.
            // It is not possible to also maintain the odd function (hence the unspecified sign).
            return result.apply(real, Double.NaN);
        }
        if (imaginary == 0) {
            // Real-only sinh(x).
            return result.apply(Math.sinh(real), imaginary);
        }
        final double x = Math.abs(real);
        if (x > SAFE_EXP) {
            // Approximate sinh/cosh(x) using exp^|x| / 2
            return coshsinh(x, real, imaginary, true, result);
        }
        // No overflow of sinh/cosh
        return result.apply(Math.sinh(real) * Math.cos(imaginary),
            Math.cosh(real) * Math.sin(imaginary));
    }
    /**
     * Returns the hyperbolic cosine of the complex number.
     *
     * <p>This function exists to allow implementation of the identity
     * {@code cos(z) = cosh(iz)}.<p>
     *
     * @param real Real part.
     * @param imaginary Imaginary part.
     * @param result Constructor.
     * @return The hyperbolic cosine of the complex number.
     */
    private static <R> R cosh(double real, double imaginary, ComplexResult<R> result) {
        // ISO C99: Preserve the even function by mapping to positive
        // f(z) = f(-z)
        if (Double.isInfinite(real) && !Double.isFinite(imaginary)) {
            return result.apply(Math.abs(real), Double.NaN);
        }
        if (real == 0) {
            // Imaginary-only cosh(iy) = cos(y).
            if (Double.isFinite(imaginary)) {
                // Maintain periodic property with respect to the imaginary component.
                // sinh(+/-0.0) * sin(+/-x) = +/-0 * sin(x)
                return result.apply(Math.cos(imaginary),
                    changeSign(real, Math.sin(imaginary)));
            }
            // If imaginary is inf/NaN the sign of the imaginary part is unspecified.
            // Although not required by C99 changing the sign maintains the conjugate equality.
            // It is not possible to also maintain the even function (hence the unspecified sign).
            return result.apply(Double.NaN, changeSign(real, imaginary));
        }
        if (imaginary == 0) {
            // Real-only cosh(x).
            // Change sign to preserve conjugate equality and even function.
            // sin(+/-0) * sinh(+/-x) = +/-0 * +/-a (sinh is monotonic and same sign)
            // => change the sign of imaginary using real. Handles special case of infinite real.
            // If real is NaN the sign of the imaginary part is unspecified.
            return result.apply(Math.cosh(real), changeSign(imaginary, real));
        }
        final double x = Math.abs(real);
        if (x > SAFE_EXP) {
            // Approximate sinh/cosh(x) using exp^|x| / 2
            return coshsinh(x, real, imaginary, false, result);
        }
        // No overflow of sinh/cosh
        return result.apply(Math.cosh(real) * Math.cos(imaginary),
            Math.sinh(real) * Math.sin(imaginary));
    }

    /**
     * Compute cosh or sinh when the absolute real component |x| is large. In this case
     * cosh(x) and sinh(x) can be approximated by exp(|x|) / 2:
     *
     * <pre>
     * cosh(x+iy) real = (e^|x| / 2) * cos(y)
     * cosh(x+iy) imag = (e^|x| / 2) * sin(y) * sign(x)
     * sinh(x+iy) real = (e^|x| / 2) * cos(y) * sign(x)
     * sinh(x+iy) imag = (e^|x| / 2) * sin(y)
     * </pre>
     *
     * @param x Absolute real component |x|.
     * @param real Real part (x).
     * @param imaginary Imaginary part (y).
     * @param sinh Set to true to compute sinh, otherwise cosh.
     * @param result Constructor.
     * @return The hyperbolic sine/cosine of the complex number.
     */
    private static <R> R coshsinh(double x, double real, double imaginary, boolean sinh,
                                    ComplexResult<R> result) {
        // Always require the cos and sin.
        double re = Math.cos(imaginary);
        double im = Math.sin(imaginary);
        // Compute the correct function
        if (sinh) {
            re = changeSign(re, real);
        } else {
            im = changeSign(im, real);
        }
        // Multiply by (e^|x| / 2).
        // Overflow safe computation since sin/cos can be very small allowing a result
        // when e^x overflows: e^x / 2 = (e^m / 2) * e^m * e^(x-2m)
        if (x > SAFE_EXP * 3) {
            // e^x > e^m * e^m * e^m
            // y * (e^m / 2) * e^m * e^m will overflow when starting with Double.MIN_VALUE.
            // Note: Do not multiply by +inf to safeguard against sin(y)=0.0 which
            // will create 0 * inf = nan.
            re *= Double.MAX_VALUE * Double.MAX_VALUE * Double.MAX_VALUE;
            im *= Double.MAX_VALUE * Double.MAX_VALUE * Double.MAX_VALUE;
        } else {
            // Initial part of (e^x / 2) using (e^m / 2)
            re *= EXP_M / 2;
            im *= EXP_M / 2;
            double xm;
            if (x > SAFE_EXP * 2) {
                // e^x = e^m * e^m * e^(x-2m)
                re *= EXP_M;
                im *= EXP_M;
                xm = x - SAFE_EXP * 2;
            } else {
                // e^x = e^m * e^(x-m)
                xm = x - SAFE_EXP;
            }
            final double exp = Math.exp(xm);
            re *= exp;
            im *= exp;
        }
        return result.apply(re, im);
    }
    /**
     * Compute {@code x^2 + y^2 - 1} in high precision.
     * Assumes that the values x and y can be multiplied without overflow; that
     * {@code x >= y}; and both values are positive.
     *
     * @param x the x value
     * @param y the y value
     * @return {@code x^2 + y^2 - 1}.
     */
    private static double x2y2m1(double x, double y) {
        // Hull et al used (x-1)*(x+1)+y*y.
        // From the paper on page 236:

        // If x == 1 there is no cancellation.

        // If x > 1, there is also no cancellation, but the argument is now accurate
        // only to within a factor of 1 + 3 EPSILSON (note that x – 1 is exact),
        // so that error = 3 EPSILON.

        // If x < 1, there can be serious cancellation:

        // If 4 y^2 < |x^2 – 1| the cancellation is not serious ... the argument is accurate
        // only to within a factor of 1 + 4 EPSILSON so that error = 4 EPSILON.

        // Otherwise there can be serious cancellation and the relative error in the real part
        // could be enormous.

        final double xx = x * x;
        final double yy = y * y;
        // Modify to use high precision before the threshold set by Hull et al.
        // This is to preserve the monotonic output of the computation at the switch.
        // Set the threshold when x^2 + y^2 is above 0.5 thus subtracting 1 results in a number
        // that can be expressed with a higher precision than any number in the range 0.5-1.0
        // due to the variable exponent used below 0.5.
        if (x < 1 && xx + yy > 0.5) {
            // Large relative error.
            // This does not use o.a.c.numbers.LinearCombination.value(x, x, y, y, 1, -1).
            // It is optimised knowing that:
            // - the products are squares
            // - the final term is -1 (which does not require split multiplication and addition)
            // - The answer will not be NaN as the terms are not NaN components
            // - The order is known to be 1 > |x| >= |y|
            // The squares are computed using a split multiply algorithm and
            // the summation using an extended precision summation algorithm.

            // Split x and y as one 26 bits number and one 27 bits number
            final double xHigh = splitHigh(x);
            final double xLow  = x - xHigh;
            final double yHigh = splitHigh(y);
            final double yLow  = y - yHigh;

            // Accurate split multiplication x * x and y * y
            final double x2Low = squareLow(xLow, xHigh, xx);
            final double y2Low = squareLow(yLow, yHigh, yy);

            return sumx2y2m1(xx, x2Low, yy, y2Low);
        }
        return (x - 1) * (x + 1) + yy;
    }

    /**
     * Implement Dekker's method to split a value into two parts. Multiplying by (2^s + 1) create
     * a big value from which to derive the two split parts.
     * <pre>
     * c = (2^s + 1) * a
     * a_big = c - a
     * a_hi = c - a_big
     * a_lo = a - a_hi
     * a = a_hi + a_lo
     * </pre>
     *
     * <p>The multiplicand must be odd allowing a p-bit value to be split into
     * (p-s)-bit value {@code a_hi} and a non-overlapping (s-1)-bit value {@code a_lo}.
     * Combined they have (p􏰔-1) bits of significand but the sign bit of {@code a_lo}
     * contains a bit of information.
     *
     * @param a Value.
     * @return the high part of the value.
     * @see <a href="https://doi.org/10.1007/BF01397083">
     * Dekker (1971) A floating-point technique for extending the available precision</a>
     */
    private static double splitHigh(double a) {
        final double c = MULTIPLIER * a;
        return c - (c - a);
    }

    /**
     * Compute the round-off from the square of a split number with {@code low} and {@code high}
     * components. Uses Dekker's algorithm for split multiplication modified for a square product.
     *
     * <p>Note: This is candidate to be replaced with {@code Math.fma(x, x, -x * x)} to compute
     * the round-off from the square product {@code x * x}. This would remove the requirement
     * to compute the split number and make this method redundant. {@code Math.fma} requires
     * JDK 9 and FMA hardware support.
     *
     * @param low Low part of number.
     * @param high High part of number.
     * @param square Square of the number.
     * @return <code>low * low - (((product - high * high) - low * high) - high * low)</code>
     * @see <a href="http://www-2.cs.cmu.edu/afs/cs/project/quake/public/papers/robust-arithmetic.ps">
     * Shewchuk (1997) Theorum 18</a>
     */
    private static double squareLow(double low, double high, double square) {
        final double lh = low * high;
        return low * low - (((square - high * high) - lh) - lh);
    }

    /**
     * Compute the round-off from the sum of two numbers {@code a} and {@code b} using
     * Dekker's two-sum algorithm. The values are required to be ordered by magnitude:
     * {@code |a| >= |b|}.
     *
     * @param a First part of sum.
     * @param b Second part of sum.
     * @param x Sum.
     * @return <code>b - (x - a)</code>
     * @see <a href="http://www-2.cs.cmu.edu/afs/cs/project/quake/public/papers/robust-arithmetic.ps">
     * Shewchuk (1997) Theorum 6</a>
     */
    private static double fastSumLow(double a, double b, double x) {
        // x = a + b
        // bVirtual = x - a
        // y = b - bVirtual
        return b - (x - a);
    }

    /**
     * Compute the round-off from the sum of two numbers {@code a} and {@code b} using
     * Knuth's two-sum algorithm. The values are not required to be ordered by magnitude.
     *
     * @param a First part of sum.
     * @param b Second part of sum.
     * @param x Sum.
     * @return <code>(a - (x - (x - a))) + (b - (x - a))</code>
     * @see <a href="http://www-2.cs.cmu.edu/afs/cs/project/quake/public/papers/robust-arithmetic.ps">
     * Shewchuk (1997) Theorum 7</a>
     */
    private static double sumLow(double a, double b, double x) {
        // x = a + b
        // bVirtual = x - a
        // aVirtual = x - bVirtual
        // bRoundoff = b - bVirtual
        // aRoundoff = a - aVirtual
        // y = aRoundoff + bRoundoff
        final double bVirtual = x - a;
        return (a - (x - bVirtual)) + (b - bVirtual);
    }

    /**
     * Sum x^2 + y^2 - 1. It is assumed that {@code y <= x < 1}.
     *
     * <p>Implement Shewchuk's expansion-sum algorithm: [x2Low, x2High] + [-1] + [y2Low, y2High].
     *
     * @param x2High High part of x^2.
     * @param x2Low Low part of x^2.
     * @param y2High High part of y^2.
     * @param y2Low Low part of y^2.
     * @return x^2 + y^2 - 1
     * @see <a href="http://www-2.cs.cmu.edu/afs/cs/project/quake/public/papers/robust-arithmetic.ps">
     * Shewchuk (1997) Theorum 12</a>
     */
    private static double sumx2y2m1(double x2High, double x2Low, double y2High, double y2Low) {
        // Let e and f be non-overlapping expansions of components of length m and n.
        // The following algorithm will produce a non-overlapping expansion h where the
        // sum h_i = e + f and components of h are in increasing order of magnitude.

        // Expansion-sum proceeds by a grow-expansion of the first part from one expansion
        // into the other, extending its length by 1. The process repeats for the next part
        // but the grow-expansion starts at the previous merge position + 1.
        // Thus expansion-sum requires mn two-sum operations to merge length m into length n
        // resulting in length m+n-1.

        // Variables numbered from 1 as per Figure 7 (p.12). The output expansion h is placed
        // into e increasing its length for each grow expansion.

        // We have two expansions for x^2 and y^2 and the whole number -1.
        // Expecting (x^2 + y^2) close to 1 we generate first the intermediate expansion
        // (x^2 - 1) moving the result away from 1 where there are sparse floating point
        // representations. This is then added to a similar magnitude y^2. Leaving the -1
        // until last suffers from 1 ulp rounding errors more often and the requirement
        // for a distillation sum to reduce rounding error frequency.

        // Note: Do not use the alternative fast-expansion-sum of the parts sorted by magnitude.
        // The parts can be ordered with a single comparison into:
        // [y2Low, (y2High|x2Low), x2High, -1]
        // The fast-two-sum saves 1 fast-two-sum and 3 two-sum operations (21 additions) and
        // adds a penalty of a single branch condition.
        // However the order in not "strongly non-overlapping" and the fast-expansion-sum
        // output will not be strongly non-overlapping. The sum of the output has 1 ulp error
        // on random cis numbers approximately 1 in 160 events. This can be removed by a
        // distillation two-sum pass over the final expansion as a cost of 1 fast-two-sum and
        // 3 two-sum operations! So we use the expansion sum with the same operations and
        // no branches.

        // q=running sum
        double q = x2Low - 1;
        double e1 = fastSumLow(-1, x2Low, q);
        double e3 = q + x2High;
        double e2 = sumLow(q, x2High, e3);

        final double f1 = y2Low;
        final double f2 = y2High;

        // Grow expansion of f1 into e
        q = f1 + e1;
        e1 = sumLow(f1, e1, q);
        double p = q + e2;
        e2 = sumLow(q, e2, p);
        double e4 = p + e3;
        e3 = sumLow(p, e3, e4);

        // Grow expansion of f2 into e (only required to start at e2)
        q = f2 + e2;
        e2 = sumLow(f2, e2, q);
        p = q + e3;
        e3 = sumLow(q, e3, p);
        final double e5 = p + e4;
        e4 = sumLow(p, e4, e5);

        // Final summation:
        // The sum of the parts is within 1 ulp of the true expansion value e:
        // |e - sum| < ulp(sum).
        // To achieve the exact result requires iteration of a distillation two-sum through
        // the expansion until convergence, i.e. no smaller term changes higher terms.
        // This requires (n-1) iterations for length n. Here we neglect this as
        // although the method is not ensured to be exact is it robust on random
        // cis numbers.
        return e1 + e2 + e3 + e4 + e5;
    }
}
