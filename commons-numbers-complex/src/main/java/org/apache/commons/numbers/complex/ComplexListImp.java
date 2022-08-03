package org.apache.commons.numbers.complex;

import java.util.Arrays;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.ListIterator;

public class ComplexListImp implements List<ComplexDouble> {
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
    protected double[] realAndImagParts;
    /** TODO. */
    private int size;

    public ComplexListImp() {
        this(DEFAULT_CAPACITY);
    }

    public ComplexListImp(int capacity) {
        capacity = Math.max(DEFAULT_CAPACITY, capacity);
        realAndImagParts = new double[capacity];
    }

    protected ComplexListImp(double[] realAndImag, int size) {
        realAndImagParts = realAndImag;
        this.size = size;
    }

    /** TODO. */
    @Override
    public int size() {
        return size;
    }

    @Override
    public boolean isEmpty() {
        return size == 0;
    }

    @Override
    public boolean contains(Object o) {
        return false;
    }

    @Override
    public Iterator<ComplexDouble> iterator() {
        return null;
    }

    @Override
    public Object[] toArray() {
        return new Object[0];
    }

    @Override
    public <T> T[] toArray(T[] a) {
        return null;
    }


    /** TODO.
     * @param startIndex
     * @param length
     * @param copy
     * @return dest
     **/
    protected double[] getDestinationRealPart(int startIndex, int length, boolean copy) {
        return this.realAndImagParts;
    }

    /** TODO.
     * @param startIndex
     * @param length
     * @param copy
     * @return dest
     **/
    protected double[] getDestinationImaginaryPart(int startIndex, int length, boolean copy) {
        return this.realAndImagParts;
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
    protected  ComplexListImp getComplexListImp(double[] real, double[] img, int length) {
        return this;
    }

    private void grow() {
        grow(size + 1);
    }

    private void grow(int minCapacity) {
        // overflow-conscious code
        int oldCapacity = realAndImagParts.length;
        int newCapacity = oldCapacity + (oldCapacity >> 1);
        if (newCapacity - minCapacity < 0) {
            newCapacity = minCapacity;
        }
        if (newCapacity - MAX_ARRAY_SIZE > 0) {
            newCapacity = hugeCapacity(minCapacity);
        }
        // minCapacity is usually close to size, so this is a win:
        realAndImagParts = Arrays.copyOf(realAndImagParts, newCapacity);
        realAndImagParts = Arrays.copyOf(realAndImagParts, newCapacity);
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


    @Override
    public final ComplexDouble get(int index) {
        return Complex.ofCartesian(realAndImagParts[index], realAndImagParts[index+1]);
    }
    

    @Override
    public final ComplexDouble set(int index, ComplexDouble element) {
        rangeCheck(index);

        Complex oldValue = Complex.ofCartesian(realAndImagParts[index], realAndImagParts[index]);
        realAndImagParts[index] = element.real();
        realAndImagParts[index] = element.imag();
        return oldValue;

    }

    @Override
    public final ComplexDouble remove(int index) {
        rangeCheck(index);

        modCount++;
        Complex oldValue = Complex.ofCartesian(realAndImagParts[index], realAndImagParts[index]);

        int numMoved = size - index - 1;
        if (numMoved > 0) {
            System.arraycopy(realAndImagParts, index + 1, realAndImagParts, index,
                numMoved);
            System.arraycopy(realAndImagParts, index + 1, realAndImagParts, index,
                numMoved);
        }
        // clear to let GC do its work
        realAndImagParts[--size] = 0;
        realAndImagParts[--size] = 0;

        return oldValue;
    }

    @Override
    public int indexOf(Object o) {
        return 0;
    }

    @Override
    public int lastIndexOf(Object o) {
        return 0;
    }

    @Override
    public ListIterator<ComplexDouble> listIterator() {
        return null;
    }

    @Override
    public ListIterator<ComplexDouble> listIterator(int index) {
        return null;
    }

    @Override
    public List<ComplexDouble> subList(int fromIndex, int toIndex) {
        return null;
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

    @Override
    public boolean remove(Object o) {
        return false;
    }

    @Override
    public boolean containsAll(Collection<?> c) {
        return false;
    }

    @Override
    public boolean addAll(Collection<? extends ComplexDouble> c) {
        return false;
    }

    @Override
    public boolean addAll(int index, Collection<? extends ComplexDouble> c) {
        return false;
    }

    @Override
    public boolean removeAll(Collection<?> c) {
        return false;
    }

    @Override
    public boolean retainAll(Collection<?> c) {
        return false;
    }

    public final boolean add(double real, double imag) {
        modCount++;
        final int s;
        if ((s = size) == (this.realAndImagParts).length) {
            grow();
        }
        this.realAndImagParts[size] = real;
        this.realAndImagParts[size] = imag;
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
        if ((s = size) == (this.realAndImagParts).length) {
            grow();
        }
        System.arraycopy(this.realAndImagParts, index,
            this.realAndImagParts, index + 1,
            s - index);
        this.realAndImagParts[index] = element.real();
        this.realAndImagParts[index+1] = element.imag();
        size = s + 1;
    }



    @Override
    public final boolean equals(Object other) {
        if (this == other) {
            return true;
        }

        if (other instanceof ComplexListImp) {
            final ComplexListImp c = (ComplexListImp) other;
            if (this.size != c.size) {
                return false;
            }
            return equals(this.realAndImagParts, c.realAndImagParts, size);
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
        result = 31 * result + hashCode(this.realAndImagParts, size);
        return 31 * result + hashCode(this.realAndImagParts, size);
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
                    .append(this.realAndImagParts[i]).append(FORMAT_SEP)
                    .append(this.realAndImagParts[i])
                    .append(FORMAT_END)
                    .append(ARRAY_SEP);
            }
            builder.append(FORMAT_START)
                .append(this.realAndImagParts[len]).append(FORMAT_SEP)
                .append(this.realAndImagParts[len])
                .append(FORMAT_END);
        }
        return builder.append(ARRAY_FORMAT_END).toString();
    }

    public final String toString(int index) {
        rangeCheck(index);
        return new StringBuilder(TO_STRING_SIZE)
            .append(FORMAT_START)
            .append(this.realAndImagParts[index]).append(FORMAT_SEP)
            .append(this.realAndImagParts[index])
            .append(FORMAT_END).toString();
    }


    public final Complex forEach(int index, ComplexUnaryOperator operator) {
        return (Complex) forEach(index, operator, ComplexSink.D_COMPLEX_RESULT);
    }

    public final ComplexDouble forEach(int index, ComplexUnaryOperator operator, ComplexSink<ComplexDouble> result) {
        this.rangeCheck(index);
        return operator.apply(this.realAndImagParts[index], this.realAndImagParts[index], result);
    }

    public final Complex forEach(int index, Complex operand, ComplexBinaryOperator operator) {
        return (Complex) forEach(index, operand, operator, ComplexSink.D_COMPLEX_RESULT);
    }

    public final ComplexDouble forEach(int index, Complex operand, ComplexBinaryOperator operator,
                                       ComplexSink<ComplexDouble> result) {
        this.rangeCheck(index);
        return operator.apply(this.realAndImagParts[index], this.realAndImagParts[index], operand.real(),
            operand.imag(), result);
    }


    public final ComplexListImp forEach(ComplexUnaryOperator operator) {
        return  forEach(0, size, operator);
    }
    public final ComplexListImp forEach(int startIndex, int length, ComplexUnaryOperator operator) {
        rangeCheckForSubList(startIndex, length);
        double[] destinationRealPart = getDestinationRealPart(startIndex, length, false);
        double[] destinationImaginaryPart = getDestinationImaginaryPart(startIndex, length, false);
        int s = getDestinationStartIndex(startIndex, length);
        ComplexListImp.Cursor cursor = new ComplexListImp.Cursor(destinationRealPart, destinationImaginaryPart,
            realAndImagParts, realAndImagParts, startIndex, s);
        for (int i = 0; i < length; i++) {
            cursor.setIndex(i);
            operator.apply(cursor.real(), cursor.imag(), cursor);
        }
        return getComplexListImp(destinationRealPart, destinationImaginaryPart, length);
    }

    public final ComplexListImp forEach(Complex operand, ComplexBinaryOperator operator) {
        return  forEach(operand.real(), operand.imag(), operator);
    }

    public final ComplexListImp forEach(double realOperand, double imaginaryOperand,
                                     ComplexBinaryOperator operator) {
        return  forEach(0, size, realOperand,  imaginaryOperand, operator);
    }

    public final ComplexListImp forEach(int startIndex, int length, Complex operand,
                                     ComplexBinaryOperator operator) {
        return  forEach(startIndex, length, operand.real(), operand.imag(),
            operator);
    }

    private static final class Cursor  implements ComplexDouble, ComplexSink<Void> {
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
    public final ComplexListImp forEach(int startIndex, int length, double realOperand,
                                     double imaginaryOperand, ComplexBinaryOperator operator) {
        rangeCheckForSubList(startIndex, length);
        double[] destinationRealPart = getDestinationRealPart(startIndex, length, false);
        double[] destinationImaginaryPart = getDestinationImaginaryPart(startIndex, length, false);
        int s = getDestinationStartIndex(startIndex, length);
        ComplexListImp.Cursor cursor = new ComplexListImp.Cursor(destinationRealPart, destinationImaginaryPart,
            realAndImagParts, realAndImagParts, startIndex, s);

        for (int i = 0; i < length; i++) {
            cursor.setIndex(i);
            operator.apply(cursor.real(), cursor.imag(), realOperand,
                imaginaryOperand, cursor);
        }
        return getComplexListImp(destinationRealPart, destinationImaginaryPart, length);
    }
    
}
