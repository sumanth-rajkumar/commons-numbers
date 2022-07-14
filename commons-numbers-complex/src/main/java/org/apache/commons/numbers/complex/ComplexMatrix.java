package org.apache.commons.numbers.complex;

import java.io.Serializable;

public interface ComplexMatrix extends Serializable {

    /**
     * Returns the number of rows in this matrix.
     *
     * @return Number of rows.
     */
    int getNumRows();

    /**
     * Returns the number of columns in this matrix.
     *
     * @return Number of columns.
     */
    int getNumCols();

    /**
     * Returns the complex value of the matrix's element
     * @param row Matrix element's row index..
     * @param col Matrix element's column index.
     */
    ComplexDouble get(int row, int col);

    /**
     * Set's the complex value of the matrix's element
     *
     * @param row Matrix element's row index..
     * @param col Matrix element's column index.
     * @param real The real component
     * @param imaginary The imaginary component
     */
    void set(int row, int col, double real, double imaginary);

    /**
     * Returns the real component of the matrix's element.
     *
     * @param row Matrix element's row index..
     * @param col Matrix element's column index.
     * @return The specified element's value.
     */
    double getReal(int row, int col);


    /**
     * Sets the real component of the matrix's element.
     *
     * @param row Matrix element's row index..
     * @param col Matrix element's column index.
     * @param val  The element's new value.
     */
    void setReal(int row, int col, double val);

    /**
     * Returns the imaginary component of the matrix's element.
     *
     * @param row Matrix element's row index..
     * @param col Matrix element's column index.
     * @return The specified element's value.
     */
    double getImag(int row, int col);


    /**
     * Sets the imaginary component of the matrix's element.
     *
     * @param row Matrix element's row index..
     * @param col Matrix element's column index.
     * @param val  The element's new value.
     */
    void setImag(int row, int col, double val);

    /**
     * Returns the number of elements in the internal data array
     *
     * @return Number of elements in the data array.
     */
    int getDataLength();

    ComplexVector getRow(int row);

    ComplexVector getCol(int col);

}
