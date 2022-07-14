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
import java.util.Iterator;
import java.util.List;



public class ComplexList extends AbstractList<ComplexDouble> implements List<ComplexDouble>, ComplexVector {
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

    @Override
    public void get(int index, int destIndex, int len, double[] realAndImgPairs) {

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
    public final ComplexDouble get(int index) {
        return Complex.ofCartesian(realParts[index], imaginaryParts[index]);
    }

    @Override
    public ComplexVector set(int index, double real, double imaginary) {
        realParts[index] = real;
        imaginaryParts[index] = imaginary;
        return this;
    }

    @Override
    public double getReal(int i) {
        return realParts[i];
    }

    @Override
    public ComplexVector setReal(int i, double val) {
        return null;
    }

    @Override
    public double getImag(int i) {
        return imaginaryParts[i];
    }

    @Override
    public ComplexVector setImag(int i, double val) {
        return null;
    }

    /**
     * To do.
     * @param index
     * @param length
     * @return
    */
    @Override
    public Iterator<ComplexDouble> iterator(int index, int length) {
        //To do.
        throw new UnsupportedOperationException();
    }

    /**
     * To do.
     * @param index
     * @param sourceIndex
     * @param len
     * @param realAndImgPairs
     * @return
    */
    @Override
    public ComplexVector setValues(int index, int sourceIndex, int len, double[] realAndImgPairs) {
        for (int i = 0; i < len; i += 2) {
            final int srcOffset = sourceIndex + i * 2;
            this.realParts[index + i] = realAndImgPairs[srcOffset];
            this.imaginaryParts[index + i] = realAndImgPairs[srcOffset + 1];
        }
        return this;
    }

    /**
     * To do.
     * @param index
     * @param c
     * @return
    */
    @Override
    public ComplexVector setValue(int index, ComplexDouble c) {
        realParts[index] = c.real();
        imaginaryParts[index] = c.imag();
        return this;
    }

    /**
     * To do.
     * @param index
     * @param r
     * @param i
     * @return
    */
    @Override
    public ComplexVector setValue(int index, double r, double i) {
        realParts[index] = r;
        imaginaryParts[index] = i;
        return this;
    }

    @Override
    public final ComplexDouble set(int index, ComplexDouble element) {
        rangeCheck(index);

        Complex oldValue = Complex.ofCartesian(realParts[index], imaginaryParts[index]);
        realParts[index] = element.real();
        imaginaryParts[index] = element.imag();
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
    public final boolean add(ComplexDouble element) {
        return add(element.real(), element.imag());
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
    public final void add(int index, ComplexDouble element) {
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
        this.realParts[index] = element.real();
        this.imaginaryParts[index] = element.imag();
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

    // public double real()
    public final double real(int index) {
        rangeCheck(index);
        return realParts[index];
    }

    public final double[] realList(int index, int length) {
        rangeCheckForSubList(index, length);
        double[] result = new double[length];
        System.arraycopy(realParts, index, result, 0, length);
        return result;
    }

    public final double imag(int index) {
        rangeCheck(index);
        return imaginaryParts[index];
    }

    public final double[] imagList(int index, int length) {
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


    public final Complex forEach(int index, ComplexUnaryOperator operator) {
        return (Complex) forEach(index, operator, ComplexConstructor.D_COMPLEX_RESULT);
    }

    public final ComplexDouble forEach(int index, ComplexUnaryOperator operator, ComplexConstructor<ComplexDouble> result) {
        this.rangeCheck(index);
        return operator.apply(this.realParts[index], this.imaginaryParts[index], result);
    }

    public final Complex forEach(int index, Complex operand, ComplexBinaryOperator operator) {
        return (Complex) forEach(index, operand, operator, ComplexConstructor.D_COMPLEX_RESULT);
    }

    public final ComplexDouble forEach(int index, Complex operand, ComplexBinaryOperator operator,
                                       ComplexConstructor<ComplexDouble> result) {
        this.rangeCheck(index);
        return operator.apply(this.realParts[index], this.imaginaryParts[index], operand.real(),
            operand.imag(), result);
    }


    public final ComplexList forEach(ComplexUnaryOperator operator) {
        return  forEach(0, size, operator);
    }
    public final ComplexList forEach(int startIndex, int length, ComplexUnaryOperator operator) {
        rangeCheckForSubList(startIndex, length);
        double[] destinationRealPart = getDestinationRealPart(startIndex, length, false);
        double[] destinationImaginaryPart = getDestinationImaginaryPart(startIndex, length, false);
        int s = getDestinationStartIndex(startIndex, length);
        Cursor cursor = new Cursor(destinationRealPart, destinationImaginaryPart,
            realParts, imaginaryParts, startIndex, s);
        for (int i = 0; i < length; i++) {
            cursor.setIndex(i);
            operator.apply(cursor.real(), cursor.imag(), cursor);
        }
        return getComplexList(destinationRealPart, destinationImaginaryPart, length);
    }

    public final ComplexList forEach(Complex operand, ComplexBinaryOperator operator) {
        return  forEach(operand.real(), operand.imag(), operator);
    }

    public final ComplexList forEach(double realOperand, double imaginaryOperand,
                                     ComplexBinaryOperator operator) {
        return  forEach(0, size, realOperand,  imaginaryOperand, operator);
    }

    public final ComplexList forEach(int startIndex, int length, Complex operand,
                                     ComplexBinaryOperator operator) {
        return  forEach(startIndex, length, operand.real(), operand.imag(),
            operator);
    }

    private static final class Cursor  implements ComplexDouble, ComplexConstructor<Void> {
        /**
         * To do.
        */
        private int offset;

        /**
         * To do.
         */
        private int destStart;

        /**
         * To do.
         */
        private int srcStart;
        /**
         * To do.
        */
        private double[] destinationRealPart;

        /**
         * To do.
        */
        private double[] destinationImgPart;

        /**
         * To do.
         */
        private double[] srcRealPart;

        /**
         * To do.
         */
        private double[] srcImgPart;

        private Cursor(double[] r, double[] i, double[] sr, double[] si, int s, int d) {
            destinationRealPart = r;
            destinationImgPart = i;
            srcRealPart = sr;
            srcImgPart = si;
            srcStart = s;
            destStart = d;
        }

        @Override
        public Void apply(double r, double i) {
            int of = destStart + offset;
            destinationRealPart[of] = r;
            destinationImgPart[of] = i;
            return null;
        }

        public void setIndex(int i) {
            offset = i;
        }

        @Override
        public double real() {
            return srcRealPart[srcStart + offset];
        }

        @Override
        public double imag() {
            return srcImgPart[srcStart + offset];
        }
    }
    public final ComplexList forEach(int startIndex, int length, double realOperand,
                                     double imaginaryOperand, ComplexBinaryOperator operator) {
        rangeCheckForSubList(startIndex, length);
        double[] destinationRealPart = getDestinationRealPart(startIndex, length, false);
        double[] destinationImaginaryPart = getDestinationImaginaryPart(startIndex, length, false);
        int s = getDestinationStartIndex(startIndex, length);
        Cursor cursor = new Cursor(destinationRealPart, destinationImaginaryPart,
            realParts, imaginaryParts, startIndex, s);

        for (int i = 0; i < length; i++) {
            cursor.setIndex(i);
            operator.apply(cursor.real(), cursor.imag(), realOperand,
                                                     imaginaryOperand, cursor);
        }
        return getComplexList(destinationRealPart, destinationImaginaryPart, length);
    }

    public final Complex conj(int index) {
        return forEach(index, ComplexFunctions::conj);
    }

    public final ComplexList conj() {
        //return (ComplexList) apply(ComplexFunctions::conj);
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
        return forEach(index, factor, Complex::multiply);
    }

    public final ComplexList multiply(Complex factor) {
        return forEach(0, size, factor.real(), factor.imag(), ComplexFunctions::multiply);
    }

    public final ComplexList multiply(int startIndex, int length, Complex factor) {
        return forEach(startIndex, length, factor.real(), factor.imag(), ComplexFunctions::multiply);
    }

    public final Complex divide(int index, Complex factor) {
        return forEach(index, factor, Complex::divide);
    }

    public final ComplexList divide(Complex factor) {
        return forEach(0, size, factor.real(), factor.imag(), ComplexFunctions::divide);
    }

    public final ComplexList divide(int startIndex, int length, Complex factor) {
        return forEach(startIndex, length, factor.real(), factor.imag(), ComplexFunctions::divide);
    }


}


