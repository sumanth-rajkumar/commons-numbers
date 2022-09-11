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

package org.apache.commons.numbers.complex.arrays;

import org.apache.commons.numbers.complex.Complex;
import org.apache.commons.numbers.complex.ComplexConsumer;
import org.apache.commons.numbers.complex.ComplexUnaryOperator;

import java.util.AbstractList;
import java.util.Arrays;
import java.util.Collection;
import java.util.ConcurrentModificationException;
import java.util.Objects;

/**
 * This is an abstract class that implements all the common data and operations for a ComplexList.
 *
 * <p>The memory usage is more efficient than using a List of Complex objects as the underlying numbers are not stored
 * using instances of Complex.</p>
 *
 * <p>This list does not support {@code null} Complex objects.
 */
public abstract class ComplexList extends AbstractList<Complex> {

    /**
     * The maximum size of array to allocate.
     * Ensuring max capacity is even with additional space for VM array headers.
     */
    protected static final int MAX_ARRAY_SIZE = Integer.MAX_VALUE - 9;
    /** Default initial capacity. */
    protected static final int DEFAULT_CAPACITY = 8;
    /** Memory error message. */
    protected static final String OOM_ERROR_STRING = "cannot allocate capacity %s greater than max ";
    /** Size label message. */
    private static final String SIZE_MSG = ", Size: ";
    /** Index position label message. */
    private static final String INDEX_MSG = "Index: ";

    /** Size of complex numbers in the list. */
    protected int size;

    /**
     * Constructor to prevent extension of this class outside inner clases.
     */
    private ComplexList() {
        //Intentionally empty
    }

    /** {@inheritDoc} */
    @Override
    public int size() {
        return size;
    }

    /**
     * Checks if the given index is in range.
     *
     * @param index Index of the complex number to range check.
     * @throws IndexOutOfBoundsException if index isn't within the range.
     */
    protected void rangeCheck(int index) {
        if (index >= size) {
            throw new IndexOutOfBoundsException(outOfBoundsMsg(index));
        }
    }

    /**
     * A version of rangeCheck used by add and addAll.
     *
     * @param index Index of the complex number to range check.
     * @throws IndexOutOfBoundsException if index isn't within the range of list.
     */
    protected void rangeCheckForInsert(int index) {
        if (index < 0 || index > size) {
            throw new IndexOutOfBoundsException(outOfBoundsMsg(index));
        }
    }

    /**
     * Gets the complex number \( (a + i b) \) at the indexed position of the list.
     *
     * @return the complex number.
     */
    @Override
    public abstract Complex get(int index);

    /**
     * Gets the real part \( a \) of the complex number \( (a + i b) \) at the indexed position of the list.
     *
     * @param index Index of the complex number.
     * @return The real part.
     */
    public abstract double getReal(int index);

    /**
     * Gets the imaginary part \( b \) of the complex number \( (a + i b) \) at the indexed position of the list.
     *
     * @param index Index of the complex number.
     * @return The imaginary part.
     */
    public abstract double getImaginary(int index);

    /**
     * Replaces the complex number at the specified position in this list with the specified
     * complex number.
     *
     * @param index Index of the complex number.
     * @param number Complex number to be set.
     * @return The previous number.
     * @throws NullPointerException if the number is null.
     */
    @Override
    public abstract Complex set(int index, Complex number);

    /**
     * Replaces the complex number's real part at the specified position
     * in the list with the specified real number.
     *
     * @param index Index of the complex number.
     * @param real Real part \( a \) of the complex number \( (a +ib) \).
     */
    public abstract void setReal(int index, double real);

    /**
     * Replaces the complex number's imaginary part at the specified position
     * in the list with the specified imaginary number.
     *
     * @param index Index of the complex number.
     * @param imaginary Imaginary part \( b \) of the complex number \( (a +ib) \).
     */
    public abstract void setImaginary(int index, double imaginary);

    /**
     * Returns an array containing all the complex number's real parts in
     * the list in proper sequence (from first to last number).
     *
     * @return Array of real parts.
     */
    abstract double[] toArrayReal();

    /**
     * Returns an array containing all the complex number's imaginary parts in
     * the list in proper sequence (from first to last number).
     *
     * @return Array of imaginary parts.
     */
    abstract double[] toArrayImaginary();

    /**
     * Appends the specified complex number to the end of this list.
     *
     * @param number Complex number to be appended to this list
     * @throws NullPointerException if the number is null.
     */
    @Override
    public abstract boolean add(Complex number);

    /**
     * Inserts the specified complex number at the specified position in this list.
     *
     * @param index Index at which the specified number is to be inserted.
     * @param number Complex number to be inserted into this list.
     * @throws NullPointerException if the number is null.
     */
    @Override
    public abstract void add(int index, Complex number);

    /**
     * Appends all of the elements in the specified collection to the end of this list.
     *
     * @param c Collection containing numbers to be added to this list.
     * @throws NullPointerException if any number in the collection null.
     */
    @Override
    public abstract boolean addAll(Collection<? extends Complex> c);

    /**
     * Inserts all of the elements in the specified collection into this list at the specified position.
     *
     * @param index Index at which the specified collection is to be inserted.
     * @param c Collection containing numbers to be inserted into this list.
     * @throws NullPointerException if any number in the collection null.
     */
    @Override
    public abstract boolean addAll(int index, Collection<? extends Complex> c);

    /**
     * Removes the complex number at the specified position in this list.
     * Shifts any subsequent numbers to the left (subtracts one from their indices).
     * Returns the number that was removed from the list.
     *
     * @param index Index of the number to be removed.
     * @return The number previously at the index.
     */
    @Override
    public abstract Complex remove(int index);

    /**
     * Replaces each complex number of the list with the result of applying the operator to that complex number.
     *
     * @param operator The operator to apply to each complex number.
     * @throws NullPointerException if the specified operator is null.
     */
    public abstract void replaceAll(ComplexUnaryOperator<Void> operator);

    /**
     * Performs the given action for each complex number of the list until all complex numbers
     * have been processed or the action throws an exception. Unless otherwise specified by the
     * implementing class, actions are performed in the order of iteration.
     *
     * @param action The action to apply to each complex number.
     * @throws NullPointerException if the specified action is null.
     */
    public abstract void forEach(ComplexConsumer action);

    /**
     * Constructs an IndexOutOfBoundsException detail message.
     *
     * @param index Index of the complex number.
     * @return message detailing the exception.
     */
    private String outOfBoundsMsg(int index) {
        return INDEX_MSG + index + SIZE_MSG + size;
    }

    /**
     * Constructs an empty interleaved list which can store up to the specified capacity without a memory reallocation.
     *
     * @param capacity Capacity of interleaved list.
     * @return ComplexList object.
     * @throws IllegalArgumentException if the {@code capacity} is greater than {@code MAX_CAPACITY}.
     */
    public static ComplexList interleaved(int capacity) {
        return new ComplexInterleavedList(capacity);
    }

    /**
     * Constructs an empty non-interleaved list which can store up to the specified capacity without a memory reallocation.
     *
     * @param capacity Capacity of non-interleaved list.
     * @return ComplexList object.
     * @throws IllegalArgumentException if the {@code capacity} is greater than {@code MAX_CAPACITY}.
     */
    public static ComplexList nonInterleaved(int capacity) {
        return new ComplexNonInterleavedList(capacity);
    }

    /**
     * Constructs an empty interleaved list.
     *
     * @return ComplexList object.
     */
    public static ComplexList interleaved() {
        return new ComplexInterleavedList();
    }

    /**
     * Constructs an empty non-interleaved list.
     *
     * @return ComplexList object.
     */
    public static ComplexList nonInterleaved() {
        return new ComplexNonInterleavedList();
    }

    /**
     * Constructs an interleaved list using the specified double array.
     * The data isn't defensively copied, the specified array is used in-place
     * and therefore any external modifications to the array will reflect on the list
     * unless a structural modification like resize is made to the data storage.
     *
     * @param data Specified backing double array.
     * @return ComplexList object.
     * @throws IllegalArgumentException if the specified double array doesn't have an even amount of doubles.
     */
    public static ComplexList from(double[] data) {
        return new ComplexInterleavedList(data);
    }

    /**
     * Constructs a non-interleaved list using the specified double arrays.
     * The data isn't defensively copied, the specified arrays is used in-place
     * and therefore any external modifications to the arrays will reflect on the list
     * unless a structural modification like resize is made to the data storage.
     *
     * @param realData Specified backing double array for real parts.
     * @param imaginaryData Specified backing double array for imaginary parts.
     * @return ComplexList object.
     * @throws IllegalArgumentException if the specified double arrays don't have the same amount of doubles.
     */
    public static ComplexList from(double[] realData, double[] imaginaryData) {
        return new ComplexNonInterleavedList(realData, imaginaryData);
    }

    /**
     * Resizable-double array implementation of the List interface. Implements all optional list operations,
     * and permits all complex numbers. In addition to implementing the List interface,
     * this class provides methods to manipulate the size of the array that is used internally to store the list.
     *
     * <p>Each ComplexInterleavedList instance has a capacity. The capacity is half the size of the double array used to store the complex numbers
     * in the list. As complex numbers are added to an ComplexInterleavedList, its capacity grows automatically.
     * The complex number is stored using an interleaved format and so the maximum number of complex numbers that may be added is
     * approximately 2<sup>30</sup>. This is half the maximum capacity of java.util.ArrayList.
     * The memory usage is more efficient than using a List of Complex objects as the underlying numbers are not stored
     * using instances of Complex.</p>
     *
     * <p>An application can increase the capacity of an ComplexInterleavedList instance before adding a large number of complex numbers
     * using the ensureCapacityInternal operation. This may reduce the amount of incremental reallocation.</p>
     *
     * <p>This list does not support {@code null} Complex objects.
     */
    private static class ComplexInterleavedList extends ComplexList {

        /** Max capacity for size of complex numbers in the list. */
        private static final int MAX_CAPACITY = MAX_ARRAY_SIZE / 2;
        /** Error in case of allocation above max capacity. */
        private static final String OOM_ERROR = OOM_ERROR_STRING + MAX_CAPACITY;
        /**
         * Storage for the complex numbers.
         * Data is stored in an interleaved format using (real, imaginary) pairs.
         */
        private double[] realAndImagParts;

        /**
         * Constructs an empty list which can store up to the specified capacity without a memory reallocation.
         *
         * @param capacity Capacity of list.
         * @throws IllegalArgumentException if the {@code capacity} is greater than {@code MAX_CAPACITY}.
         */
        ComplexInterleavedList(int capacity) {
            if (capacity > MAX_CAPACITY) {
                throw new IllegalArgumentException(String.format(OOM_ERROR, capacity));
            }
            final int arrayLength = Math.max(DEFAULT_CAPACITY, capacity) * 2;
            realAndImagParts = new double[arrayLength];
        }

        /**
         * Constructs an empty list.
         */
        ComplexInterleavedList() {
            realAndImagParts = new double[DEFAULT_CAPACITY * 2];
        }

        /**
         * Constructs an interleaved list using the specified double array.
         * The data isn't defensively copied, the specified array is used in-place.
         *
         * @param fromArray Specified backing double array.
         * @throws IllegalArgumentException if the specified double array doesn't have an even amount of doubles.
         * @throws NullPointerException if the specified double array is null.
         */
        ComplexInterleavedList(double[] fromArray) {
            if ((fromArray.length & 1) != 0) {
                throw new IllegalArgumentException("Length of array has to be even");
            }
            realAndImagParts = fromArray;
            size = fromArray.length >> 1;
        }

        /** {@inheritDoc} */
        @Override
        public Complex get(int index) {
            rangeCheck(index);
            final int i = index << 1;
            return Complex.ofCartesian(realAndImagParts[i], realAndImagParts[i + 1]);
        }

        /** {@inheritDoc} */
        @Override
        public double getReal(int index) {
            rangeCheck(index);
            return realAndImagParts[index << 1];
        }

        /** {@inheritDoc} */
        @Override
        public double getImaginary(int index) {
            rangeCheck(index);
            return realAndImagParts[(index << 1) + 1];
        }

        /** {@inheritDoc} */
        @Override
        public Complex set(int index, Complex number) {
            rangeCheck(index);
            final int i = index << 1;
            final Complex oldValue = Complex.ofCartesian(realAndImagParts[i], realAndImagParts[i + 1]);
            realAndImagParts[i] = number.getReal();
            realAndImagParts[i + 1] = number.getImaginary();
            return oldValue;
        }

        /** {@inheritDoc} */
        @Override
        public void setReal(int index, double real) {
            rangeCheck(index);
            realAndImagParts[index << 1] = real;
        }

        /** {@inheritDoc} */
        @Override
        public void setImaginary(int index, double imaginary) {
            rangeCheck(index);
            realAndImagParts[(index << 1) + 1] = imaginary;
        }

        /** {@inheritDoc} */
        @Override
        double[] toArrayReal() {
            final int length = size;
            double[] real = new double[length];
            for (int i = 0; i < length; i++) {
                real[i] = realAndImagParts[i << 1];
            }
            return real;
        }

        /** {@inheritDoc} */
        @Override
        double[] toArrayImaginary() {
            final int length = size;
            double[] imag = new double[length];
            for (int i = 0; i < length; i++) {
                imag[i] = realAndImagParts[(i << 1) + 1];
            }
            return imag;
        }

        /**
         * Increases the capacity of this ComplexInterleavedList instance, if necessary, to ensure that it can hold at
         * least the amount of complex numbers specified by the minimum capacity argument.
         *
         * @param minCapacity Desired minimum capacity.
         * @return the backing double array.
         * @throws OutOfMemoryError if the {@code minCapacity} is greater than {@code MAX_ARRAY_SIZE}.
         */
        private double[] ensureCapacityInternal(int minCapacity) {
            modCount++;
            final long minArrayCapacity = Integer.toUnsignedLong(minCapacity) << 1;
            if (minArrayCapacity > MAX_ARRAY_SIZE) {
                throw new OutOfMemoryError(String.format(OOM_ERROR, minArrayCapacity));
            }
            final long oldArrayCapacity = realAndImagParts.length;
            if (minArrayCapacity > oldArrayCapacity) {
                long newArrayCapacity = oldArrayCapacity + (oldArrayCapacity >> 1);
                // Round-odd up to even
                newArrayCapacity += newArrayCapacity & 1;

                // Ensure minArrayCapacity <= newArrayCapacity <= MAX_ARRAY_SIZE
                // Note: At this point minArrayCapacity <= MAX_ARRAY_SIZE
                if (newArrayCapacity > MAX_ARRAY_SIZE) {
                    newArrayCapacity = MAX_ARRAY_SIZE;
                } else if (newArrayCapacity < minArrayCapacity) {
                    newArrayCapacity = minArrayCapacity;
                }
                realAndImagParts = Arrays.copyOf(realAndImagParts, (int) newArrayCapacity);
            }
            return realAndImagParts;
        }

        /**
         * Increases the capacity of this ComplexInterleavedList instance, if necessary, to ensure that it can hold at
         * least an additional amount of complex numbers specified by the capacity argument.
         *
         * @param capacity Desired capacity.
         */
        private void expand(int capacity) {
            ensureCapacityInternal(size + capacity);
        }

        /** {@inheritDoc} */
        @Override
        public boolean add(Complex number) {
            double[] e = realAndImagParts;
            if (size == (e.length >>> 1)) {
                e = ensureCapacityInternal(size + 1);
            }
            final int i = size << 1;
            e[i] = number.real();
            e[i + 1] = number.imag();
            size++;
            return true;
        }

        /** {@inheritDoc} */
        @Override
        public void add(int index, Complex number) {
            rangeCheckForInsert(index);
            double[] e = realAndImagParts;
            if (size == e.length >>> 1) {
                e = ensureCapacityInternal(size + 1);
            }
            final double real = number.real();
            final double imaginary = number.imag();
            final int i = index << 1;
            final int s = size << 1;
            System.arraycopy(e, i, e, i + 2, s - i);
            e[i] = real;
            e[i + 1] = imaginary;
            size++;
        }

        /** {@inheritDoc} */
        @Override
        public boolean addAll(Collection<? extends Complex> c) {
            final int numNew = c.size();
            expand(numNew);
            double[] realAndImgData = new double[numNew * 2];
            int i = 0;
            for (final Complex val : c) {
                realAndImgData[i++] = val.getReal();
                realAndImgData[i++] = val.getImaginary();
            }
            final int s = size << 1;
            System.arraycopy(realAndImgData, 0, realAndImagParts, s, realAndImgData.length);
            size += numNew;
            return numNew != 0;
        }

        /** {@inheritDoc} */
        @Override
        public boolean addAll(int index, Collection<? extends Complex> c) {
            rangeCheckForInsert(index);
            final int numNew = c.size();
            final int numNew2 = numNew << 1;
            expand(numNew);
            final double[] realAndImgData = new double[numNew * 2];
            int i = 0;
            for (final Complex val : c) {
                realAndImgData[i++] = val.getReal();
                realAndImgData[i++] = val.getImaginary();
            }
            final int numMoved = (size - index) * 2;
            final int index2 = index << 1;
            System.arraycopy(realAndImagParts, index2, realAndImagParts, index2 + numNew2, numMoved);
            System.arraycopy(realAndImgData, 0, realAndImagParts, index2, realAndImgData.length);
            size += numNew;
            return numNew != 0;
        }

        /** {@inheritDoc} */
        @Override
        public Complex remove(int index) {
            rangeCheck(index);
            modCount++;
            final int i = index << 1;
            final int s = size << 1;
            final Complex oldValue = Complex.ofCartesian(realAndImagParts[i], realAndImagParts[i + 1]);
            final int numMoved = s - i - 2;
            if (numMoved > 0) {
                System.arraycopy(realAndImagParts, i + 2, realAndImagParts, i, numMoved);
            }
            size--;
            return oldValue;
        }

        /** {@inheritDoc} */
        @Override
        public void replaceAll(ComplexUnaryOperator<Void> operator) {
            Objects.requireNonNull(operator);
            final double[] parts = this.realAndImagParts;
            final int m = size;
            final int expectedModCount = modCount;
            for (int i = 0; i < m; i++) {
                final int index = i << 1;
                operator.apply(parts[index], parts[index + 1], (x, y) -> {
                    parts[index] = x;
                    parts[index + 1] = y;
                    return null;
                });
            }
            // check for comodification
            if (modCount != expectedModCount) {
                throw new ConcurrentModificationException();
            }
            modCount++;
        }

        /** {@inheritDoc} */
        @Override
        public void forEach(ComplexConsumer action) {
            Objects.requireNonNull(action);
            final double[] parts = this.realAndImagParts;
            final int m = size;
            for (int i = 0; i < m; i++) {
                final int index = i << 1;
                action.accept(parts[index], parts[index + 1]);
            }
        }
    }

    /**
     * Resizable-double arrays implementation of the List interface. Implements all optional list operations,
     * and permits all complex numbers. In addition to implementing the List interface,
     * this class provides methods to manipulate the size of the arrays that are used internally to store the list.
     *
     * <p>Each ComplexNonInterleavedList instance has a capacity. The capacity is the size of the double arrays used to store the complex numbers
     * in the list. As complex numbers are added to an ComplexNonInterleavedList, its capacity grows automatically.
     * The complex number is stored using an non-interleaved format and so the maximum number of complex numbers that may be added is
     * approximately 2<sup>31</sup>. This is also the maximum capacity of java.util.ArrayList.
     * The memory usage is more efficient than using a List of Complex objects as the underlying numbers are not stored
     * using instances of Complex.</p>
     *
     * <p>An application can increase the capacity of an ComplexNonInterleavedList instance before adding a large number of complex numbers
     * using the ensureCapacityInternal operation. This may reduce the amount of incremental reallocation.</p>
     *
     * <p>This list does not support {@code null} Complex objects.
     */
    private static class ComplexNonInterleavedList extends ComplexList {

        /** Max capacity for size of complex numbers in the list. */
        private static final int MAX_CAPACITY = MAX_ARRAY_SIZE;
        /** Error in case of allocation above max capacity. */
        private static final String OOM_ERROR = OOM_ERROR_STRING + MAX_CAPACITY;
        /**
         * Storage for the real parts of complex numbers.
         */
        private double[] realParts;

        /**
         * Storage for the imaginary parts of complex numbers.
         */
        private double[] imaginaryParts;

        /**
         * Constructs an empty list which can store up to the specified capacity without a memory reallocation.
         *
         * @param capacity Capacity of list.
         * @throws IllegalArgumentException if the {@code capacity} is greater than {@code MAX_CAPACITY}.
         */
        ComplexNonInterleavedList(int capacity) {
            if (capacity > MAX_CAPACITY) {
                throw new IllegalArgumentException(String.format(OOM_ERROR, capacity));
            }
            final int arrayLength = Math.max(DEFAULT_CAPACITY, capacity);
            realParts = new double[arrayLength];
            imaginaryParts = new double[arrayLength];
        }

        /**
         * Constructs an empty list.
         */
        ComplexNonInterleavedList() {
            realParts = new double[DEFAULT_CAPACITY];
            imaginaryParts = new double[DEFAULT_CAPACITY];
        }

        /**
         * Constructs a non-interleaved list using the specified double arrays.
         * The data isn't defensively copied, the specified arrays is used in-place.
         *
         * @param fromReal Specified backing double array for real parts.
         * @param fromImaginary Specified backing double array for imaginary parts.
         * @throws IllegalArgumentException if the specified double arrays don't have the same amount of doubles.
         * @throws NullPointerException if either of the specified double arrays are null.
         */
        ComplexNonInterleavedList(double[] fromReal, double[] fromImaginary) {
            if (fromReal.length != fromImaginary.length) {
                throw new IllegalArgumentException("Need the same amount of real and imaginary parts");
            }
            realParts = fromReal;
            imaginaryParts = fromImaginary;
            size = fromReal.length;
        }

        /** {@inheritDoc} */
        @Override
        public Complex get(int index) {
            rangeCheck(index);
            return Complex.ofCartesian(realParts[index], imaginaryParts[index]);
        }

        /** {@inheritDoc} */
        @Override
        public double getReal(int index) {
            rangeCheck(index);
            return realParts[index];
        }

        /** {@inheritDoc} */
        @Override
        public double getImaginary(int index) {
            rangeCheck(index);
            return imaginaryParts[index];
        }

        /** {@inheritDoc} */
        @Override
        public Complex set(int index, Complex number) {
            rangeCheck(index);
            final Complex oldValue = Complex.ofCartesian(realParts[index], imaginaryParts[index]);
            realParts[index] = number.getReal();
            imaginaryParts[index] = number.getImaginary();
            return oldValue;
        }

        /** {@inheritDoc} */
        @Override
        public void setReal(int index, double real) {
            rangeCheck(index);
            realParts[index] = real;
        }

        /** {@inheritDoc} */
        @Override
        public void setImaginary(int index, double imaginary) {
            rangeCheck(index);
            imaginaryParts[index] = imaginary;
        }

        /** {@inheritDoc} */
        @Override
        double[] toArrayReal() {
            return Arrays.copyOf(realParts, size);
        }

        /** {@inheritDoc} */
        @Override
        double[] toArrayImaginary() {
            return Arrays.copyOf(imaginaryParts, size);
        }

        /**
         * Increases the capacity of this ComplexNonInterleavedList instance, if necessary, to ensure that it can hold at
         * least the amount of complex numbers specified by the minimum capacity argument.
         *
         * @param minCapacity Desired minimum capacity.
         * @throws OutOfMemoryError if the {@code minCapacity} is greater than {@code MAX_ARRAY_SIZE}.
         */
        private void ensureCapacityInternal(int minCapacity) {
            modCount++;
            final long minArrayCapacity = Integer.toUnsignedLong(minCapacity);
            if (minArrayCapacity > MAX_ARRAY_SIZE) {
                throw new OutOfMemoryError(String.format(OOM_ERROR, minArrayCapacity));
            }
            final long oldArrayCapacity = realParts.length;
            if (minArrayCapacity > oldArrayCapacity) {
                long newArrayCapacity = oldArrayCapacity + (oldArrayCapacity >> 1);

                // Ensure minArrayCapacity <= newArrayCapacity <= MAX_ARRAY_SIZE
                // Note: At this point minArrayCapacity <= MAX_ARRAY_SIZE
                if (newArrayCapacity > MAX_ARRAY_SIZE) {
                    newArrayCapacity = MAX_ARRAY_SIZE;
                } else if (newArrayCapacity < minArrayCapacity) {
                    newArrayCapacity = minArrayCapacity;
                }
                realParts = Arrays.copyOf(realParts, (int) newArrayCapacity);
                imaginaryParts = Arrays.copyOf(imaginaryParts, (int) newArrayCapacity);
            }
        }

        /**
         * Increases the capacity of this ComplexNonInterleavedList instance, if necessary, to ensure that it can hold at
         * least an additional amount of complex numbers specified by the capacity argument.
         *
         * @param capacity Desired capacity.
         */
        private void expand(int capacity) {
            ensureCapacityInternal(size + capacity);
        }

        /** {@inheritDoc} */
        @Override
        public boolean add(Complex number) {
            if (size == realParts.length) {
                ensureCapacityInternal(size + 1);
            }
            final int i = size;
            realParts[i] = number.real();
            imaginaryParts[i] = number.imag();
            size++;
            return true;
        }

        /** {@inheritDoc} */
        @Override
        public void add(int index, Complex number) {
            rangeCheckForInsert(index);
            if (size == realParts.length) {
                ensureCapacityInternal(size + 1);
            }
            final double real = number.real();
            final double imaginary = number.imag();
            final int s = size;
            System.arraycopy(realParts, index, realParts, index + 1, s - index);
            System.arraycopy(imaginaryParts, index, imaginaryParts, index + 1, s - index);
            realParts[index] = real;
            imaginaryParts[index] = imaginary;
            size++;
        }

        /** {@inheritDoc} */
        @Override
        public boolean addAll(Collection<? extends Complex> c) {
            final int numNew = c.size();
            expand(numNew);
            double[] realData = new double[numNew];
            double[] imaginaryData = new double[numNew];
            int i = 0;
            for (final Complex val : c) {
                realData[i] = val.getReal();
                imaginaryData[i] = val.getImaginary();
                i++;
            }
            final int s = size;
            System.arraycopy(realData, 0, realParts, s, realData.length);
            System.arraycopy(imaginaryData, 0, imaginaryParts, s, imaginaryData.length);
            size += numNew;
            return numNew != 0;
        }

        /** {@inheritDoc} */
        @Override
        public boolean addAll(int index, Collection<? extends Complex> c) {
            rangeCheckForInsert(index);
            final int numNew = c.size();
            expand(numNew);
            double[] realData = new double[numNew];
            double[] imaginaryData = new double[numNew];
            int i = 0;
            for (final Complex val : c) {
                realData[i] = val.getReal();
                imaginaryData[i] = val.getImaginary();
                i++;
            }
            final int numMoved = size - index;
            System.arraycopy(realData, index, realParts, index + numNew, numMoved);
            System.arraycopy(realParts, 0, realParts, index, realData.length);
            System.arraycopy(imaginaryData, index, imaginaryParts, index + numNew, numMoved);
            System.arraycopy(imaginaryParts, 0, imaginaryParts, index, imaginaryData.length);
            size += numNew;
            return numNew != 0;
        }

        /** {@inheritDoc} */
        @Override
        public Complex remove(int index) {
            rangeCheck(index);
            modCount++;
            final int s = size;
            final Complex oldValue = Complex.ofCartesian(realParts[index], imaginaryParts[index]);
            final int numMoved = s - index - 1;
            if (numMoved > 0) {
                System.arraycopy(realParts, index + 1, realParts, index, numMoved);
                System.arraycopy(imaginaryParts, index + 1, imaginaryParts, index, numMoved);
            }
            size--;
            return oldValue;
        }

        /** {@inheritDoc} */
        @Override
        public void replaceAll(ComplexUnaryOperator<Void> operator) {
            Objects.requireNonNull(operator);
            final double[] realData = this.realParts;
            final double[] imaginaryData = this.imaginaryParts;
            final int m = size;
            final int expectedModCount = modCount;
            for (int i = 0; i < m; i++) {
                final int index = i;
                operator.apply(realData[i], imaginaryData[i], (x, y) -> {
                    realData[index] = x;
                    imaginaryData[index] = y;
                    return null;
                });
            }
            // check for comodification
            if (modCount != expectedModCount) {
                throw new ConcurrentModificationException();
            }
            modCount++;
        }

        /** {@inheritDoc} */
        @Override
        public void forEach(ComplexConsumer action) {
            Objects.requireNonNull(action);
            final double[] realData = this.realParts;
            final double[] imaginaryData = this.imaginaryParts;
            final int m = size;
            for (int i = 0; i < m; i++) {
                action.accept(realData[i], imaginaryData[i]);
            }
        }
    }
}
