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

import static org.apache.commons.numbers.complex.Complex.CARTESIAN_RESULT;
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
        if (index < 0 || length < 0 || index + length > size) {
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
    public final Complex set(int index, Complex element) {
        rangeCheck(index);

        Complex oldValue = Complex.ofCartesian(realParts[index], imaginaryParts[index]);
        realParts[index] = element.getReal();
        imaginaryParts[index] = element.getImaginary();
        return oldValue;

    }

    @Override
    public final Complex remove(int index) {
        rangeCheck(index);

        modCount++;
        Complex oldValue = Complex.ofCartesian(realParts[index], imaginaryParts[index]);

        int numMoved = size - index - 1;
        if (numMoved > 0) {
            System.arraycopy(realParts, index + 1, realParts, index,
                numMoved);
            System.arraycopy(imaginaryParts, index + 1, imaginaryParts, index,
                numMoved);
        }
        // clear to let GC do its work
        realParts[--size] = 0;
        imaginaryParts[--size] = 0;

        return oldValue;
    }

    @Override
    public final void clear() {
        modCount++;
        size = 0;
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


    public final Complex forEach(int index, ComplexFunction operator) {
        return forEach(index, operator, Complex.CARTESIAN_RESULT);
    }

    public final <R> R forEach(int index, ComplexFunction operator, ComplexResult<R> result) {
        this.rangeCheck(index);
        return operator.apply(this.realParts[index], this.imaginaryParts[index], result);
    }

    public final Complex forEach(int index, Complex operand, ComplexBiFunction operator) {
        return forEach(index, operand, operator, CARTESIAN_RESULT);
    }

    public final <R> R forEach(int index, Complex operand, ComplexBiFunction operator,
                               ComplexResult<R> result) {
        this.rangeCheck(index);
        return operator.apply(this.realParts[index], this.imaginaryParts[index], operand.getReal(),
            operand.getImaginary(), result);
    }


    public final ComplexList forEach(ComplexFunction operator) {
        return  forEach(0, size, operator);
    }
    public final ComplexList forEach(int startIndex, int length, ComplexFunction operator) {
        rangeCheckForSubList(startIndex, length);
        double[] destinationRealPart = getDestinationRealPart(startIndex, length, false);
        double[] destinationImaginaryPart = getDestinationImaginaryPart(startIndex, length, false);
        int s = getDestinationStartIndex(startIndex, length);
        int offset = startIndex - s;
        int len = length + s;
        for (int i = s; i < len; i++) {
            final int destinationIndex = i;
            operator.apply(this.realParts[i + offset], this.imaginaryParts[i + offset], (x, y) -> {
                destinationRealPart[destinationIndex] = x;
                destinationImaginaryPart[destinationIndex] = y;
                return null;
            });
        }
        return getComplexList(destinationRealPart, destinationImaginaryPart, length);
    }

    public final ComplexList forEach(Complex operand, ComplexBiFunction operator) {
        return  forEach(operand.getReal(), operand.getImaginary(), operator);
    }

    public final ComplexList forEach(double realOperand, double imaginaryOperand,
                                     ComplexBiFunction operator) {
        return  forEach(0, size, realOperand,  imaginaryOperand, operator);
    }

    public final ComplexList forEach(int startIndex, int length, Complex operand,
                                     ComplexBiFunction operator) {
        return  forEach(startIndex, length, operand.getReal(), operand.getImaginary(),
            operator);
    }

    public final ComplexList forEach(int startIndex, int length, double realOperand,
                                     double imaginaryOperand, ComplexBiFunction operator) {
        rangeCheckForSubList(startIndex, length);
        double[] destinationRealPart = getDestinationRealPart(startIndex, length, false);
        double[] destinationImaginaryPart = getDestinationImaginaryPart(startIndex, length, false);
        int s = getDestinationStartIndex(startIndex, length);
        int offset = startIndex - s;
        int len = length + s;
        for (int i = s; i < len; i++) {
            final int destinationIndex = i;
            operator.apply(this.realParts[i + offset], this.imaginaryParts[i + offset], realOperand,
                                                     imaginaryOperand, (x, y) -> {
                    destinationRealPart[destinationIndex] = x;
                    destinationImaginaryPart[destinationIndex] = y;
                    return null;
                });
        }
        return getComplexList(destinationRealPart, destinationImaginaryPart, length);
    }

    public final Complex conj(int index) {
        return forEach(index, ComplexFunctions::conj);
    }

    public final ComplexList conj() {
        return forEach(0, size, ComplexFunctions::conj);
    }

    public final ComplexList conj(int startIndex, int length) {
        return forEach(startIndex, length, ComplexFunctions::conj);
    }

    public final Complex exp(int index) {
        return forEach(index, ComplexFunctions::exp);
    }

    public final ComplexList exp() {
        return forEach(0, size, ComplexFunctions::exp);
    }

    public final ComplexList exp(int startIndex, int length) {
        return forEach(startIndex, length, ComplexFunctions::exp);
    }

    public final Complex asin(int index) {
        return forEach(index, ComplexFunctions::asin);
    }

    public final ComplexList asin() {
        return forEach(0, size, ComplexFunctions::asin);
    }

    public final ComplexList asin(int startIndex, int length) {
        return forEach(startIndex, length, ComplexFunctions::asin);
    }

    public final Complex multiply(int index, Complex factor) {
        return forEach(index, factor, ComplexFunctions::multiply);
    }

    public final ComplexList multiply(Complex factor) {
        return forEach(0, size, factor.getReal(), factor.getImaginary(), ComplexFunctions::multiply);
    }

    public final ComplexList multiply(int startIndex, int length, Complex factor) {
        return forEach(startIndex, length, factor.getReal(), factor.getImaginary(), ComplexFunctions::multiply);
    }

    public final Complex divide(int index, Complex factor) {
        return forEach(index, factor, ComplexFunctions::divide);
    }

    public final ComplexList divide(Complex factor) {
        return forEach(0, size, factor.getReal(), factor.getImaginary(), ComplexFunctions::divide);
    }

    public final ComplexList divide(int startIndex, int length, Complex factor) {
        return forEach(startIndex, length, factor.getReal(), factor.getImaginary(), ComplexFunctions::divide);
    }


}


