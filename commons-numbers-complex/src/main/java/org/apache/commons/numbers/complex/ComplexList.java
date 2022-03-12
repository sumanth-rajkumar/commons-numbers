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

import java.util.AbstractList;
import java.util.Arrays;
import java.util.List;

import static org.apache.commons.numbers.complex.Complex.hypot;
import static org.apache.commons.numbers.complex.Complex.negative;

public class ComplexList extends AbstractList<Complex> implements List<Complex> {
    /** TODO. */
    private static final double[] DEFAULT_EMPTY = {};
    /** TODO. */
    private static final int DEFAULT_CAPACITY = 8;
    /** TODO. */
    private static final int MAX_ARRAY_SIZE = Integer.MAX_VALUE - 8;
    /** TODO. */
    private static final String SIZE_MSG = ", Size: ";
    /** TODO. */
    private static final String INDEX_MSG = "Index: ";
    /** TODO. **/
    private static final String ARRAY_FORMAT_START = "[";
    /** TODO. */
    private static final String ARRAY_SEP = ";";
    /** TODO. */
    private static final String ARRAY_FORMAT_END = "]";
    /** TODO. */
    private static final int TO_STRING_SIZE = 64;
    /** TODO. */
    private static final char FORMAT_START = '(';
    /** TODO. */
    private static final char FORMAT_END = ')';
    /** TODO. */
    private static final char FORMAT_SEP = ',';
    /** TODO. */
    private static final String PARSER_REGEX = "[\\[;\\]]";
    /** TODO. */
    protected double[] realParts;
    /** TODO. */
    protected double[] imaginaryParts;
    /** TODO. */
    private int size;

    public ComplexList() {
        this(DEFAULT_CAPACITY);
    }

    public ComplexList(int capacity) {
        capacity = Math.max(DEFAULT_CAPACITY, capacity);
        realParts = new double[capacity];
        imaginaryParts = new double[capacity];
    }

    protected ComplexList(double[] real, double[] img, int size) {
        realParts = real;
        imaginaryParts = img;
        this.size = size;
    }

    /** TODO. */
    @Override
    public int size() {
        return size;
    }

    /** TODO.
     * @param startIndex
     * @param length
     * @param copy
     * @return dest
     **/
    protected double[] getDestinationRealPart(int startIndex, int length, boolean copy) {
        return this.realParts;
    }

    /** TODO.
     * @param startIndex
     * @param length
     * @param copy
     * @return dest
     **/
    protected double[] getDestinationImaginaryPart(int startIndex, int length, boolean copy) {
        return this.imaginaryParts;
    }

    /** TODO.
     * @param startIndex
     * @param length
     * @return dest
     **/
    protected int getDestinationStartIndex(int startIndex, int length) {
        return startIndex;
    }

    /** TODO.
     * @param real
     * @param img
     * @param length
     * @return dest
     **/
    protected  ComplexList getComplexList(double[] real, double[] img, int length) {
        return this;
    }

    public final void ensureCapacity(int minCapacity) {
        if (minCapacity > realParts.length &&
            !(realParts == DEFAULT_EMPTY &&
            minCapacity <= DEFAULT_CAPACITY)) {
            modCount++;
            grow(minCapacity);
        }
    }

    private void grow() {
        grow(size + 1);
    }

    private void grow(int minCapacity) {
        // overflow-conscious code
        int oldCapacity = realParts.length;
        int newCapacity = oldCapacity + (oldCapacity >> 1);
        if (newCapacity - minCapacity < 0) {
            newCapacity = minCapacity;
        }
        if (newCapacity - MAX_ARRAY_SIZE > 0) {
            newCapacity = hugeCapacity(minCapacity);
        }
        // minCapacity is usually close to size, so this is a win:
        realParts = Arrays.copyOf(realParts, newCapacity);
        imaginaryParts = Arrays.copyOf(imaginaryParts, newCapacity);
    }

    private static int hugeCapacity(int minCapacity) {
        if (minCapacity < 0) {
            throw new OutOfMemoryError();
        }
        return (minCapacity > MAX_ARRAY_SIZE) ?
            Integer.MAX_VALUE :
            MAX_ARRAY_SIZE;
    }

    private String outOfBoundsMsg(int index) {
        return INDEX_MSG + index + SIZE_MSG + size;
    }

    private String outOfBoundsSubListMsg(int index, int len) {
        return INDEX_MSG + index + ", Length: " + len + SIZE_MSG + size;
    }

    private void rangeCheckForSubList(int index, int length) {
        if (index < 0 || length < 0 || index + length >= size) {
            throw new IndexOutOfBoundsException(outOfBoundsSubListMsg(index, length));
        }
    }

    private void rangeCheck(int index) {
        if (index > size || index < 0) {
            throw new IndexOutOfBoundsException(outOfBoundsMsg(index));
        }
    }

    private static void rangeCheckForRealAndImaginary(int realLen, int imgLen) {
        if (realLen != imgLen) {
            throw new RuntimeException("mismatched size realLen: " + realLen + ", imaginaryLen: " + imgLen);
        }
    }


    @Override
    public final Complex get(int index) {
        return Complex.ofCartesian(realParts[index], imaginaryParts[index]);
    }

    @Override
    public final boolean add(Complex element) {
        return add(element.getReal(), element.getImaginary());
    }

    public final boolean add(double real, double imag) {
        modCount++;
        final int s;
        if ((s = size) == (this.realParts).length) {
            grow();
        }
        this.realParts[size] = real;
        this.imaginaryParts[size] = imag;
        size = s + 1;
        return true;
    }

    public final boolean addPolar(double rho, double theta) {
        return add(rho * Math.cos(theta), rho * Math.sin(theta));
    }

    public final boolean addCis(double theta) {
        return add(Math.cos(theta), Math.sin(theta));
    }

    @Override
    public final void add(int index, Complex element) {
        rangeCheck(index);
        modCount++;
        final int s;
        if ((s = size) == (this.realParts).length) {
            grow();
        }
        System.arraycopy(this.realParts, index,
            this.realParts, index + 1,
            s - index);
        System.arraycopy(this.imaginaryParts, index,
            this.imaginaryParts, index + 1,
            s - index);
        this.realParts[index] = element.getReal();
        this.imaginaryParts[index] = element.getImaginary();
        size = s + 1;
    }


    public final void addAll(double[] real, double[] imaginary) {
        rangeCheckForRealAndImaginary(real.length, imaginary.length);

        if (real.length == 0) {
            return;
        }

        modCount++;
        final int s = size;
        if (this.realParts.length - size < real.length) {
            grow(real.length + size);
        }
        System.arraycopy(real, 0, this.realParts, size, real.length);
        System.arraycopy(imaginary, 0, this.imaginaryParts, size, imaginary.length);
        size += real.length;
    }

    public static ComplexList ofCartesian(double[] real, double[] imaginary) {
        rangeCheckForRealAndImaginary(real.length, imaginary.length);

        ComplexList r = new ComplexList(real.length);
        if (real.length == 0) {
            return r;
        }

        r.addAll(real, imaginary);
        return r;
    }

    public static ComplexList ofPolar(double[] rho, double[] theta) {
        rangeCheckForRealAndImaginary(rho.length, theta.length);

        ComplexList r = new ComplexList(rho.length);
        if (rho.length == 0) {
            return r;
        }

        for (int i = 0; i < rho.length; i++) {
            // Require finite theta and non-negative, non-nan rho
            if (!Double.isFinite(theta[i]) || negative(rho[i]) || Double.isNaN(rho[i])) {
                r.realParts[i] = Double.NaN;
                r.imaginaryParts[i] = Double.NaN;
            } else {
                r.realParts[i] = rho[i] * Math.cos(theta[i]);
                r.imaginaryParts[i] = rho[i] * Math.sin(theta[i]);
            }
        }
        r.size = rho.length;
        return r;
    }

    public static ComplexList ofCis(double[] theta) {
        ComplexList r = new ComplexList(theta.length);
        if (theta.length == 0) {
            return r;
        }

        for (int i = 0; i < theta.length; i++) {
            // Require finite theta and non-negative, non-nan rho
            if (!Double.isFinite(theta[i])) {
                r.realParts[i] = Double.NaN;
                r.imaginaryParts[i] = Double.NaN;
            } else {
                r.realParts[i] =  Math.cos(theta[i]);
                r.imaginaryParts[i] =  Math.sin(theta[i]);
            }
        }
        r.size = theta.length;
        return r;
    }

    public static ComplexList parse(String s) {
        String[] strings = s.split(PARSER_REGEX);
        ComplexList r = new ComplexList(strings.length);
        for (String str : strings) {
            if (str.isEmpty()) {
                continue;
            }
            Complex com = Complex.parse(str);
            r.add(com);
        }
        return r;
    }

    // public double getReal()
    public final double getReal(int index) {
        rangeCheck(index);
        return realParts[index];
    }

    public final double[] getRealList(int index, int length) {
        rangeCheckForSubList(index, length);
        double[] result = new double[length];
        System.arraycopy(realParts, index, result, 0, length);
        return result;
    }

    public final double getImaginary(int index) {
        rangeCheck(index);
        return imaginaryParts[index];
    }

    public final double[] getImaginaryList(int index, int length) {
        rangeCheckForSubList(index, length);
        double[] result = new double[length];
        System.arraycopy(imaginaryParts, index, result, 0, length);
        return result;
    }

    public final double abs(int index) {
        rangeCheck(index);
        return hypot(this.realParts[index], this.imaginaryParts[index]);
    }

    public final double[] absList(int index, int length) {
        rangeCheckForSubList(index, length);
        double[] absResult = new double[length];
        for (int i = 0; i < length; i++) {
            absResult[i] = hypot(this.realParts[index + i], this.imaginaryParts[index + i]);
        }
        return absResult;
    }

    public final double arg(int index) {
        rangeCheck(index);
        return Math.atan2(this.imaginaryParts[index], this.realParts[index]);
    }

    public final double[] argList(int index, int length) {
        rangeCheckForSubList(index, length);
        double[] argResult = new double[length];
        for (int i = 0; i < length; i++) {
            argResult[i] = Math.atan2(this.imaginaryParts[index + i], this.realParts[index + i]);
        }
        return argResult;
    }

    public final double norm(int index) {
        rangeCheck(index);
        if (isInfinite(index)) {
            return Double.POSITIVE_INFINITY;
        }
        return this.realParts[index] * this.realParts[index] + this.imaginaryParts[index] * this.imaginaryParts[index];
    }

    public final boolean isNaN(int index) {
        rangeCheck(index);
        if (Double.isNaN(this.realParts[index]) || Double.isNaN(this.imaginaryParts[index])) {
            return !isInfinite(index);
        }
        return false;
    }

    public final boolean isInfinite(int index) {
        rangeCheck(index);
        return Double.isInfinite(this.realParts[index]) || Double.isInfinite(this.imaginaryParts[index]);
    }

    public final boolean isFinite(int index) {
        rangeCheck(index);
        return Double.isFinite(this.realParts[index]) && Double.isFinite(this.imaginaryParts[index]);
    }

    public final ComplexList conj() {
        return this.conj(0, size);
    }

    public final Complex conj(int index) {
        rangeCheck(index);
        return Complex.ofCartesian(realParts[index], -imaginaryParts[index]);
    }

    public final ComplexList conj(int startIndex, int length) {
        rangeCheckForSubList(startIndex, length);
        double[] destinationRealPart = getDestinationRealPart(startIndex, length, true);
        double[] destinationImaginaryPart = getDestinationRealPart(startIndex, length, false);

        int s = getDestinationStartIndex(startIndex, length);
        int offset = startIndex - s;
        int len = length + s;
        for (int i = s; i < len; i++) {
            destinationImaginaryPart[i] = -this.imaginaryParts[i + offset];
        }
        return getComplexList(destinationRealPart, destinationImaginaryPart, length);
    }

    public final ComplexList negate() {
        return this.negate(0, size);
    }

    public final Complex negate(int index) {
        rangeCheck(index);
        return Complex.ofCartesian(-realParts[index], -imaginaryParts[index]);
    }

    public final ComplexList negate(int startIndex, int length) {
        rangeCheckForSubList(startIndex, length);
        double[] destinationRealPart = getDestinationRealPart(startIndex, length, false);
        double[] destinationImaginaryPart = getDestinationRealPart(startIndex, length, false);
        int s = getDestinationStartIndex(startIndex, length);
        int offset = startIndex - s;
        int len = length + s;
        for (int i = s; i < len; i++) {
            destinationRealPart[i] = -this.realParts[i + offset];
            destinationImaginaryPart[i] = -this.imaginaryParts[i + offset];
        }
        return getComplexList(destinationRealPart, destinationImaginaryPart, length);
    }
    public final ComplexList proj() {
        return this.proj(0, size);
    }

    public final Complex proj(int index) {
        rangeCheck(index);
        if (isInfinite(index)) {
            return Complex.ofCartesian(Double.POSITIVE_INFINITY, Math.copySign(0.0, imaginaryParts[index]));
        }
        return Complex.ofCartesian(realParts[index], imaginaryParts[index]);
    }

    public final ComplexList proj(int startIndex, int length) {
        rangeCheckForSubList(startIndex, length);
        double[] destinationRealPart = getDestinationRealPart(startIndex, length, false);
        double[] destinationImaginaryPart = getDestinationRealPart(startIndex, length, false);
        int s = getDestinationStartIndex(startIndex, length);
        int offset = startIndex - s;
        int len = length + s;
        for (int i = s; i < len; i++) {
            if (isInfinite(i + offset)) {
                destinationRealPart[i] = Double.POSITIVE_INFINITY;
                destinationImaginaryPart[i] = Math.copySign(0.0,   this.imaginaryParts[i + offset]);
            } else {
                destinationRealPart[i] = this.realParts[i + offset];
                destinationImaginaryPart[i] = this.imaginaryParts[i + offset];
            }
        }
        return getComplexList(destinationRealPart, destinationImaginaryPart, length);
    }

    public final Complex addition(int index, Complex addend) {
        rangeCheck(index);
        return Complex.ofCartesian(this.realParts[index] + addend.real(),
            this.imaginaryParts[index] + addend.imag());
    }

    public final ComplexList addition(int startIndex, int length, Complex addend) {
        return addition(startIndex, length, addend.real(), addend.imag());
    }

    public final ComplexList addition(Complex addend) {
        return addition(0, size, addend.real(), addend.imag());
    }

    public final ComplexList addition(double realAddend, double imgAddend) {
        return addition(0, size, realAddend, imgAddend);
    }

    public final ComplexList addition(int startIndex, int length, double realAddend, double imgAddend) {
        rangeCheckForSubList(startIndex, length);
        double[] destinationRealPart = getDestinationRealPart(startIndex, length, false);
        double[] destinationImaginaryPart = getDestinationRealPart(startIndex, length, false);
        int s = getDestinationStartIndex(startIndex, length);
        int offset = startIndex - s;
        int len = length + s;
        for (int i = s; i < len; i++) {
            destinationRealPart[i] = this.realParts[i + offset] + realAddend;
            destinationImaginaryPart[i] = this.imaginaryParts[i + offset] + imgAddend;
        }
        return getComplexList(destinationRealPart, destinationImaginaryPart, length);
    }

    public final Complex addReal(int index, double addend) {
        rangeCheck(index);
        return Complex.ofCartesian(this.realParts[index] + addend,
            this.imaginaryParts[index]);
    }

    public final ComplexList addReal(double addend) {
        return addReal(0, size, addend);
    }

    public final ComplexList addReal(int startIndex, int length, double addend) {
        rangeCheckForSubList(startIndex, length);
        double[] destinationRealPart = getDestinationRealPart(startIndex, length, false);
        double[] destinationImaginaryPart = getDestinationRealPart(startIndex, length, true);
        int s = getDestinationStartIndex(startIndex, length);
        int offset = startIndex - s;
        int len = length + s;
        for (int i = s; i < len; i++) {
            destinationRealPart[i] = this.realParts[i + offset] + addend;
        }
        return getComplexList(destinationRealPart, destinationImaginaryPart, length);

    }

    public final Complex addImaginary(int index, double addend) {
        rangeCheck(index);
        return Complex.ofCartesian(this.realParts[index],
            this.imaginaryParts[index] + addend);
    }

    public final ComplexList addImaginary(double addend) {
        return  addImaginary(0, size, addend);
    }

    public final ComplexList addImaginary(int startIndex, int length, double addend) {
        rangeCheckForSubList(startIndex, length);
        double[] destinationRealPart = getDestinationRealPart(startIndex, length, true);
        double[] destinationImaginaryPart = getDestinationRealPart(startIndex, length, false);
        int s = getDestinationStartIndex(startIndex, length);
        int offset = startIndex - s;
        int len = length + s;
        for (int i = s; i < len; i++) {
            destinationImaginaryPart[i] = this.imaginaryParts[i + offset] + addend;
        }
        return getComplexList(destinationRealPart, destinationImaginaryPart, length);

    }

    public final Complex subtract(int index, Complex addend) {
        rangeCheck(index);
        return Complex.ofCartesian(this.realParts[index] - addend.real(),
            this.imaginaryParts[index] - addend.imag());
    }

    public final ComplexList subtract(int startIndex, int length, Complex addend) {
        return addition(startIndex, length, -addend.real(), -addend.imag());
    }

    public final ComplexList subtract(Complex addend) {
        return addition(0, size, -addend.real(), -addend.imag());
    }

    public final ComplexList subtract(double realAddend, double imgAddend) {
        return addition(0, size, -realAddend, -imgAddend);
    }

    public final ComplexList subtract(int startIndex, int length, double realAddend, double imgAddend) {
        return addition(startIndex, length, -realAddend, -imgAddend);
    }

    public final Complex subtract(int index, double subtrahend) {
        rangeCheck(index);
        return Complex.ofCartesian(subtrahend - this.realParts[index],
            this.imaginaryParts[index]);
    }

    public final ComplexList subtract(double subtrahend) {
        return subtract(0, size, subtrahend);
    }

    public final ComplexList subtract(int startIndex, int length, double subtrahend) {
        rangeCheckForSubList(startIndex, length);
        double[] destinationRealPart = getDestinationRealPart(startIndex, length, false);
        double[] destinationImaginaryPart = getDestinationRealPart(startIndex, length, true);
        int s = getDestinationStartIndex(startIndex, length);
        int offset = startIndex - s;
        int len = length + s;
        for (int i = s; i < len; i++) {
            destinationRealPart[i] = subtrahend - this.realParts[i + offset];
        }
        return getComplexList(destinationRealPart, destinationImaginaryPart, length);

    }

    public final Complex subtractImaginary(int index, double subtrahend) {
        rangeCheck(index);
        return Complex.ofCartesian(this.realParts[index],
            subtrahend - this.imaginaryParts[index]);
    }
    public final ComplexList subtractImaginary(double subtrahend) {
        return subtractImaginary(0, size, subtrahend);
    }

    public final ComplexList subtractImaginary(int startIndex, int length, double subtrahend) {
        rangeCheckForSubList(startIndex, length);
        double[] destinationRealPart = getDestinationRealPart(startIndex, length, true);
        double[] destinationImaginaryPart = getDestinationRealPart(startIndex, length, false);
        int s = getDestinationStartIndex(startIndex, length);
        int offset = startIndex - s;
        int len = length + s;
        for (int i = s; i < len; i++) {
            destinationImaginaryPart[i] = subtrahend - this.imaginaryParts[i + offset];
        }
        return getComplexList(destinationRealPart, destinationImaginaryPart, length);
    }

    public final Complex subtractFrom(int index, Complex subtrahend) {
        rangeCheck(index);
        return Complex.ofCartesian(subtrahend.real() - this.realParts[index],
            subtrahend.imag() - this.imaginaryParts[index]);
    }

    public final ComplexList subtractFrom(int startIndex, int length, Complex subtrahend) {
        return subtractFrom(startIndex, length, -subtrahend.real(), -subtrahend.imag());
    }

    public final ComplexList subtractFrom(Complex subtrahend) {
        return subtractFrom(0, size, -subtrahend.real(), -subtrahend.imag());
    }

    public final ComplexList subtractFrom(double realSubtrahend, double imagSubtrahend) {
        return subtractFrom(0, size, -realSubtrahend, -imagSubtrahend);
    }

    public final ComplexList subtractFrom(int startIndex, int length, double realSubtrahend, double imagSubtrahend) {
        rangeCheckForSubList(startIndex, length);
        double[] destinationRealPart = getDestinationRealPart(startIndex, length, false);
        double[] destinationImaginaryPart = getDestinationRealPart(startIndex, length, false);
        int s = getDestinationStartIndex(startIndex, length);
        int offset = startIndex - s;
        int len = length + s;
        for (int i = s; i < len; i++) {
            destinationRealPart[i] = realSubtrahend - this.realParts[i + offset];
            destinationImaginaryPart[i] = realSubtrahend - this.imaginaryParts[i + offset];
        }
        return getComplexList(destinationRealPart, destinationImaginaryPart, length);

    }

    //TODO
    public final Complex multiply(Complex factor) {
        return null;
    }

    public final Complex multiply(int index, double factor) {
        rangeCheck(index);
        return Complex.ofCartesian(this.realParts[index] * factor,
            this.imaginaryParts[index] * factor);
    }

    public final ComplexList multiply(double factor) {
        return multiply(0, size, factor);
    }

    public final ComplexList multiply(int startIndex, int length, double factor) {
        double[] destinationRealPart = getDestinationRealPart(startIndex, length, false);
        double[] destinationImaginaryPart = getDestinationRealPart(startIndex, length, false);
        int s = getDestinationStartIndex(startIndex, length);
        int offset = startIndex - s;
        int len = length + s;
        for (int i = s; i < len; i++) {
            destinationRealPart[i] = this.realParts[i + offset] * factor;
            destinationImaginaryPart[i] = this.imaginaryParts[i + offset] * factor;
        }
        return getComplexList(destinationRealPart, destinationImaginaryPart, length);

    }

    public final Complex multiplyImaginary(int index, double factor) {
        rangeCheck(index);
        return Complex.ofCartesian(-this.imaginaryParts[index] * factor, this.realParts[index] * factor);
    }

    public final ComplexList multiplyImaginary(double factor) {
        return multiplyImaginary(0, size, factor);
    }

    public final ComplexList multiplyImaginary(int startIndex, int length, double factor) {
        rangeCheckForSubList(startIndex, length);
        double[] destinationRealPart = getDestinationRealPart(startIndex, length, false);
        double[] destinationImaginaryPart = getDestinationRealPart(startIndex, length, false);
        int s = getDestinationStartIndex(startIndex, length);
        int offset = startIndex - s;
        int len = length + s;
        for (int i = s; i < len; i++) {
            destinationRealPart[i] = -this.imaginaryParts[i + offset] * factor;
            destinationImaginaryPart[i] = this.realParts[i + offset] * factor;
        }
        return getComplexList(destinationRealPart, destinationImaginaryPart, length);
    }

    //TODO
    public final Complex divide(Complex divisor) {
        return null;
    }

    public final Complex divide(int index, double divisor) {
        rangeCheck(index);
        return Complex.ofCartesian(this.realParts[index] / divisor, this.imaginaryParts[index] / divisor);
    }

    public final ComplexList divide(double divisor) {
        return divide(0, size, divisor);
    }

    public final ComplexList divide(int startIndex, int length, double divisor) {
        double[] destinationRealPart = getDestinationRealPart(startIndex, length, false);
        double[] destinationImaginaryPart = getDestinationRealPart(startIndex, length, false);
        int s = getDestinationStartIndex(startIndex, length);
        int offset = startIndex - s;
        int len = length + s;
        for (int i = s; i < len; i++) {
            destinationRealPart[i] = this.realParts[i + offset] / divisor;
            destinationImaginaryPart[i] = this.imaginaryParts[i + offset] / divisor;
        }
        return getComplexList(destinationRealPart, destinationImaginaryPart, length);

    }

    public final Complex divideImaginary(int index, double divisor) {
        rangeCheck(index);
        return Complex.ofCartesian(this.imaginaryParts[index] / divisor, -this.realParts[index] / divisor);
    }

    public final ComplexList divideImaginary(double divisor) {
        return divideImaginary(0, size, divisor);
    }

    public final ComplexList divideImaginary(int startIndex, int length, double divisor) {
        rangeCheckForSubList(startIndex, length);
        double[] destinationRealPart = getDestinationRealPart(startIndex, length, false);
        double[] destinationImaginaryPart = getDestinationRealPart(startIndex, length, false);
        int s = getDestinationStartIndex(startIndex, length);
        int offset = startIndex - s;
        int len = length + s;
        for (int i = s; i < len; i++) {
            destinationRealPart[i] = this.imaginaryParts[i + offset] / divisor;
            destinationImaginaryPart[i] = -this.realParts[i + offset] / divisor;
        }
        return getComplexList(destinationRealPart, destinationImaginaryPart, length);
    }

    //TODO
    public final Complex exp(int index) {
        return null;
    }

    public final ComplexList exp() {
        return exp(0, size);
    }

    //TODO
    public final ComplexList exp(int startIndex, int length) {
        return null;
    }

    //TODO
    public final Complex log(int index) {
        return null;
    }

    public final ComplexList log() {
        return log(0, size);
    }

    //TODO
    public final ComplexList log(int startIndex, int length) {
        return null;
    }

    //TODO
    public final Complex log10(int index) {
        return null;
    }

    public final ComplexList log10() {
        return log10(0, size);
    }

    //TODO
    public final ComplexList log10(int startIndex, int length) {
        return null;
    }

    //TODO
    public final Complex pow(int index, Complex x) {
        return null;
    }
    public final ComplexList pow(Complex x) {
        return pow(0, size, x);
    }

    //TODO
    public final ComplexList pow(int startIndex, int length, Complex x) {
        return null;
    }

    //TODO
    public final Complex pow(int index, double x) {
        return null;
    }

    public final ComplexList pow(double x) {
        return pow(0, size, x);
    }

    //TODO
    public final ComplexList pow(int startIndex, int length, double x) {
        return null;
    }

    //TODO
    public final Complex sqrt(int index) {
        return null;
    }

    public final ComplexList sqrt() {
        return sqrt(0, size);
    }

    //TODO
    public final ComplexList sqrt(int startIndex, int length) {
        return null;
    }

    //TODO
    public final Complex sin(int index) {
        return null;
    }

    public final ComplexList sin() {
        return sin(0, size);
    }

    //TODO
    public final ComplexList sin(int startIndex, int length) {
        return null;
    }

    //TODO
    public final Complex cos(int index) {
        return null;
    }

    public final ComplexList cos() {
        return cos(0, size);
    }

    //TODO
    public final ComplexList cos(int startIndex, int length) {
        return null;
    }

    //TODO
    public final Complex tan(int index) {
        return null;
    }

    public final ComplexList tan() {
        return tan(0, size);
    }

    //TODO
    public final ComplexList tan(int startIndex, int length) {
        return null;
    }

    //TODO
    public final Complex asin(int index) {
        return null;
    }

    public final ComplexList asin() {
        return asin(0, size);
    }

    //TODO
    public final ComplexList asin(int startIndex, int length) {
        return null;
    }

    //TODO
    public final Complex acos(int index) {
        return null;
    }

    public final ComplexList acos() {
        return acos(0, size);
    }

    //TODO
    public final ComplexList acos(int startIndex, int length) {
        return null;
    }

    //TODO
    public final Complex atan(int index) {
        return null;
    }

    public final ComplexList atan() {
        return atan(0, size);
    }

    //TODO
    public final ComplexList atan(int startIndex, int length) {
        return null;
    }

    //TODO
    public final Complex sinh(int index) {
        return null;
    }

    public final ComplexList sinh() {
        return sinh(0, size);
    }

    //TODO
    public final ComplexList sinh(int startIndex, int length) {
        return null;
    }

    //TODO
    public final Complex cosh(int index) {
        return null;
    }

    public final ComplexList cosh() {
        return cosh(0, size);
    }

    //TODO
    public final ComplexList cosh(int startIndex, int length) {
        return null;
    }

    //TODO
    public final Complex tanh(int index) {
        return null;
    }

    public final ComplexList tanh() {
        return tanh(0, size);
    }

    //TODO
    public final ComplexList tanh(int startIndex, int length) {
        return null;
    }

    //TODO
    public final Complex asinh(int index) {
        return null;
    }

    public final ComplexList asinh() {
        return asinh(0, size);
    }

    //TODO
    public final ComplexList asinh(int startIndex, int length) {
        return null;
    }

    //TODO
    public final Complex acosh(int index) {
        return null;
    }

    public final ComplexList acosh() {
        return acosh(0, size);
    }

    //TODO
    public final ComplexList acosh(int startIndex, int length) {
        return null;
    }

    //TODO
    public final Complex atanh(int index) {
        return null;
    }

    public final ComplexList atanh() {
        return atanh(0, size);
    }

    //TODO
    public final ComplexList atanh(int startIndex, int length) {
        return null;
    }

    //TODO
    public final ComplexList nthRoot(int index, int n) {
        return null;
    }

    @Override
    public final boolean equals(Object other) {
        if (this == other) {
            return true;
        }

        if (other instanceof ComplexList) {
            final ComplexList c = (ComplexList) other;
            if (this.size != c.size) {
                return false;
            }
            return equals(this.realParts, c.realParts, size) &&
                equals(this.imaginaryParts, c.imaginaryParts, size);
        }
        return false;
    }

    public static boolean equals(double[] a, double[] a2, int length) {
        if (a == a2) {
            return true;
        }
        if (a == null || a2 == null) {
            return false;
        }
        for (int i = 0; i < length; i++) {
            if (Double.doubleToLongBits(a[i]) != Double.doubleToLongBits(a2[i])) {
                return false;
            }
        }
        return true;
    }

    @Override
    public final int hashCode() {
        int result = 31 * size;
        result = 31 * result + hashCode(this.realParts, size);
        return 31 * result + hashCode(this.imaginaryParts, size);
    }

    public static int hashCode(double[] a, int length) {
        if (a == null) {
            return 0;
        }
        int result = 1;
        for (int i = 0; i < length; i++) {
            long bits = Double.doubleToLongBits(a[i]);
            result = 31 * result + (int)(bits ^ (bits >>> 32));
        }
        return result;
    }

    @Override
    public final String toString() {
        return toString(0, size);
    }

    public final String toString(int startIndex, int length) {
        StringBuilder builder = new StringBuilder(TO_STRING_SIZE * length);
        builder.append(ARRAY_FORMAT_START);
        if (length > 0) {
            int len = startIndex + length - 1;
            for (int i = startIndex; i < len; i++) {
                builder.append(FORMAT_START)
                    .append(this.realParts[i]).append(FORMAT_SEP)
                    .append(this.imaginaryParts[i])
                    .append(FORMAT_END)
                    .append(ARRAY_SEP);
            }
            builder.append(FORMAT_START)
                .append(this.realParts[len]).append(FORMAT_SEP)
                .append(this.imaginaryParts[len])
                .append(FORMAT_END);
        }
        return builder.append(ARRAY_FORMAT_END).toString();
    }

    public final String toString(int index) {
        rangeCheck(index);
        return new StringBuilder(TO_STRING_SIZE)
            .append(FORMAT_START)
            .append(this.realParts[index]).append(FORMAT_SEP)
            .append(this.imaginaryParts[index])
            .append(FORMAT_END).toString();
    }
}


