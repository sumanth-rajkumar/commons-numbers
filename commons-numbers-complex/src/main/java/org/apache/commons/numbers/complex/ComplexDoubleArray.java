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


public interface ComplexDoubleArray extends Iterable<ComplexDouble>, ComplexArrayResult {

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

    default ComplexDoubleArray apply(ComplexDoubleUnaryOperator op) {
        return op.apply(this, this);
    }

    default void get(double[] realAndImgPairs) {
        get(0, 0, size(), realAndImgPairs);
    }

    void get(int index, int destIndex, int len, double[] realAndImgPairs);

    ComplexDouble get(int index);

    default Iterator<ComplexDouble> iterator() {
        return iterator(0, size());
    }

    Iterator<ComplexDouble> iterator(int index, int length);

    ComplexDoubleArray setValues(int index, int sourceIndex, int len, double[] realAndImgPairs);

    ComplexDoubleArray setValue(int index, ComplexDouble c);

    ComplexDoubleArray setValue(int index, double r, double i);

    void ensureCapacity(int capacity);


}
