package org.apache.commons.numbers.complex;

import java.util.Iterator;

public interface ComplexVector {

    /**
     * Returns the complex value of the matrix's element
     * @param index vector element's index..
     */
    ComplexDouble get(int index);

    /**
     * Set's the complex value of the matrix's element
     *
     * @param i vector element's index..
     * @param real The real component
     * @param imaginary The imaginary component
     */
    ComplexVector set(int i, double real, double imaginary);

    /**
     * Returns the real component of the matrix's element.
     *
     * @param i vector element's index..
     * @return The specified element's value.
     */
    double getReal(int i);


    /**
     * Sets the real component of the matrix's element.
     *
     * @param i vector element's index..
     * @param val  The element's new value.
     */
    ComplexVector setReal(int i, double val);

    /**
     * Returns the imaginary component of the matrix's element.
     *
     * @param i vector element's index..
     * @return The specified element's value.
     */
    double getImag(int i);


    /**
     * Sets the imaginary component of the matrix's element.
     *
     * @param i vector element's index..
     * @param val  The element's new value.
     */
    ComplexVector setImag(int i, double val);

    int size();

    default boolean isImmutable() {
        return true;
    }


    default double[] toDoubleArray() {
        final int len = size();
        final double[] result = new double[len * 2];
        get(0, 0, len, result);
        return result;
    }

    default ComplexVector apply(ComplexDoubleUnaryOperator op) {
        return op.apply(this, this);
    }

    default void get(double[] realAndImgPairs) {
        get(0, 0, size(), realAndImgPairs);
    }

    void get(int index, int destIndex, int len, double[] realAndImgPairs);

    default Iterator<ComplexDouble> iterator() {
        return iterator(0, size());
    }

    Iterator<ComplexDouble> iterator(int index, int length);

    ComplexVector setValues(int index, int sourceIndex, int len, double[] realAndImgPairs);

    ComplexVector setValue(int index, ComplexDouble c);

    ComplexVector setValue(int index, double r, double i);

    void ensureCapacity(int capacity);
}
