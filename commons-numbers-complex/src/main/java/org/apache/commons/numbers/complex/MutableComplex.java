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

import java.util.Iterator;

public final class MutableComplex implements MutableComplexDouble, MutableComplexDoubleArray,
    ComplexResult<MutableComplex>, ArrayComplexResult<MutableComplex> {

    /**
     * To do.
    */
    private double real;

    /**
     * To do.
    */
    private double imaginary;

    private MutableComplex(double real, double imaginary) {
        this.real = real;
        this.imaginary = imaginary;
    }

    public static MutableComplex ofCartesian(double real, double imaginary) {
        return new MutableComplex(real, imaginary);
    }

    @Override
    public MutableComplex apply(double realP, double imaginaryP) {
        this.real = realP;
        this.imaginary = imaginaryP;
        return this;
    }

    @Override
    public MutableComplex apply(int index, double realP, double imaginaryP) {
        this.real = realP;
        this.imaginary = imaginaryP;
        return this;
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
    public void setCartesina(double realP, double imaginaryP) {
        this.real = realP;
        this.imaginary = imaginaryP;
    }

    @Override
    public int size() {
        return 1;
    }

    @Override
    public MutableComplexDoubleArray asMutable() {
        return this;
    }

    @Override
    public ComplexDoubleArray asImmutable() {
        return Complex.ofCartesian(real, imaginary);
    }

    @Override
    public void get(int index, int destIndex, int len, double[] realAndImgPairs) {
        if (index != 0 || len != 1) {
            throw new IndexOutOfBoundsException();
        }
        realAndImgPairs[destIndex] = real;
        realAndImgPairs[destIndex + 1] = imaginary;
    }

    @Override
    public MutableComplexDouble get(int index) {
        if (index != 0) {
            throw new IndexOutOfBoundsException();
        }
        return this;
    }

    @Override
    public Iterator<ComplexDouble> iterator(int index, int length) {
        return new Complex.SingletonIterator(this);
    }

    @Override
    public MutableComplexDoubleArray set(int index, int sourceIndex, int len, double[] realAndImgPairs) {
        if (index != 0 || len != 1) {
            throw new IndexOutOfBoundsException();
        }
        this.real = realAndImgPairs[sourceIndex];
        this.imaginary = realAndImgPairs[sourceIndex + 1];
        return this;
    }

    @Override
    public MutableComplexDoubleArray set(int index, ComplexDouble c) {
        if (index != 0) {
            throw new IndexOutOfBoundsException();
        }
        this.real = c.real();
        this.imaginary = c.imag();
        return this;
    }

    @Override
    public MutableComplexDoubleArray add(int index, ComplexDouble c) {
        throw new UnsupportedOperationException();
    }

    @Override
    public MutableComplexDoubleArray ensureCapacity(int capacity) {
        if (capacity != 1) {
            throw new UnsupportedOperationException();
        }
        return this;
    }
}
