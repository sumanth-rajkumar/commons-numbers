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

public class ImmutableComplexList extends ComplexList {

    protected ImmutableComplexList(double[] real, double[] imag, int size) {
        super(real, imag, size);
    }

    /** TODO.
     * @param startIndex
     * @param copy
     * @param length
     * @return dest
     **/
    @Override
    protected double[] getDestinationRealPart(int startIndex, int length, boolean copy) {
        double[] result = new double[length];
        if (copy) {
            System.arraycopy(this.realParts, startIndex, result, 0, length);
        }
        return result;
    }

    /** TODO.
     * @param startIndex
     * @param copy
     * @param length
     * @return dest
     **/
    @Override
    protected double[] getDestinationImaginaryPart(int startIndex, int length, boolean copy) {
        double[] result = new double[length];
        if (copy) {
            System.arraycopy(this.realParts, startIndex, result, 0, length);
        }
        return result;
    }

    /** TODO.
     * @param startIndex
     * @param length
     * @return dest
     **/
    @Override
    protected int getDestinationStartIndex(int startIndex, int length) {
        return 0;
    }

    /** TODO.
     * @param real
     * @param imag
     * @param length
     * @return dest
     **/
    @Override
    protected  ComplexList getComplexList(double[] real, double[] imag, int length) {
        return new ImmutableComplexList(realParts, imag, length);
    }
}
