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


public final class ComplexCartesianImpl implements Complex {
    /**
     * the real part.
     */
    private final double real;
    /**
     * the imaginary part.
     */
    private final double imaginary;

    /**
     * Private default constructor.
     *
     * @param real      Real part.
     * @param imaginary Imaginary part.
     */
    ComplexCartesianImpl(double real, double imaginary) {
        this.real = real;
        this.imaginary = imaginary;
    }

    /**
     * Test for equality with another object. If the other object is a {@code Complex} then a
     * comparison is made of the real and imaginary parts; otherwise {@code false} is returned.
     *
     * <p>If both the real and imaginary parts of two complex numbers
     * are exactly the same the two {@code Complex} objects are considered to be equal.
     * For this purpose, two {@code double} values are considered to be
     * the same if and only if the method {@link Double #doubleToLongBits(double)}
     * returns the identical {@code long} value when applied to each.
     *
     * <p>Note that in most cases, for two instances of class
     * {@code Complex}, {@code c1} and {@code c2}, the
     * value of {@code c1.equals(c2)} is {@code true} if and only if
     *
     * <pre>
     *  {@code c1.getReal() == c2.getReal() && c1.getImaginary() == c2.getImaginary()}</pre>
     *
     * <p>also has the value {@code true}. However, there are exceptions:
     *
     * <ul>
     *  <li>
     *   Instances that contain {@code NaN} values in the same part
     *   are considered to be equal for that part, even though {@code Double.NaN == Double.NaN}
     *   has the value {@code false}.
     *  </li>
     *  <li>
     *   Instances that share a {@code NaN} value in one part
     *   but have different values in the other part are <em>not</em> considered equal.
     *  </li>
     *  <li>
     *   Instances that contain different representations of zero in the same part
     *   are <em>not</em> considered to be equal for that part, even though {@code -0.0 == 0.0}
     *   has the value {@code true}.
     *  </li>
     * </ul>
     *
     * <p>The behavior is the same as if the components of the two complex numbers were passed
     * to {@link java.util.Arrays#equals(double[], double[]) Arrays.equals(double[], double[])}:
     *
     * <pre>
     *  Arrays.equals(new double[]{c1.getReal(), c1.getImaginary()},
     *                new double[]{c2.getReal(), c2.getImaginary()}); </pre>
     *
     * @param other Object to test for equality with this instance.
     * @return {@code true} if the objects are equal, {@code false} if object
     * is {@code null}, not an instance of {@code Complex}, or not equal to
     * this instance.
     * @see java.lang.Double#doubleToLongBits(double)
     * @see java.util.Arrays#equals(double[], double[])
     */
    @Override
    public  boolean equals(Object other) {
        if (this == other) {
            return true;
        }
        if (other instanceof Complex) {
            final Complex c = (Complex) other;
            return Complex.equals(real, c.real()) &&
                Complex.equals(imaginary, c.imag());
        }
        return false;
    }

    /**
     * Gets a hash code for the complex number.
     *
     * <p>The behavior is the same as if the components of the complex number were passed
     * to {@link java.util.Arrays#hashCode(double[]) Arrays.hashCode(double[])}:
     *
     * <pre>
     *  {@code Arrays.hashCode(new double[] {getReal(), getImaginary()})}</pre>
     *
     * @return A hash code value for this object.
     * @see java.util.Arrays#hashCode(double[]) Arrays.hashCode(double[])
     */
    @Override
    public int hashCode() {
        return 31 * (31 + Double.hashCode(real)) + Double.hashCode(imaginary);
    }

    /**
     * Returns a string representation of the complex number.
     *
     * <p>The string will represent the numeric values of the real and imaginary parts.
     * The values are split by a separator and surrounded by parentheses.
     * The string can be {@link #parse(String) parsed} to obtain an instance with the same value.
     *
     * <p>The format for complex number \( x + i y \) is {@code "(x,y)"}, with \( x \) and
     * \( y \) converted as if using {@link Double#toString(double)}.
     *
     * @return A string representation of the complex number.
     * @see #parse(String)
     * @see Double#toString(double)
     */
    @Override
    public String toString() {
        return new StringBuilder(TO_STRING_SIZE)
            .append(FORMAT_START)
            .append(real).append(FORMAT_SEP)
            .append(imaginary)
            .append(FORMAT_END)
            .toString();
    }

    @Override
    public double getReal() {
        return real;
    }

    @Override
    public double real() {
        return real;
    }


    @Override
    public double imag() {
        return imaginary;
    }

    @Override
    public double getImaginary() {
        return imaginary;
    }
}
