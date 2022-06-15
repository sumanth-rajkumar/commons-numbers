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

import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.simple.RandomSource;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

public class ComplexListTest {


    private ComplexList createList() {
        final double[] parts = {Double.NEGATIVE_INFINITY, -1, -0.0, 0.0, 1, Math.PI, Double.POSITIVE_INFINITY,
            Double.NaN};
        ComplexList list = new ComplexList();
        for (final double x : parts) {
            for (final double y : parts) {
                //final Complex z = Complex.ofCartesian(x, y);
                //Assertions.assertEquals(z, Complex.parse(z.toString()));
                list.add(x, y);
            }
        }
        final UniformRandomProvider rng = RandomSource.SPLIT_MIX_64.create();
        for (int i = 0; i < 10; i++) {
            final double x = -1 + rng.nextDouble() * 2;
            final double y = -1 + rng.nextDouble() * 2;
            //final Complex z = Complex.ofCartesian(x, y);
            //Assertions.assertEquals(z, Complex.parse(z.toString()));
            list.add(x, y);
        }

        // Special values not covered
        //Assertions.assertEquals(Complex.ofPolar(2, pi), Complex.parse(Complex.ofPolar(2, pi).toString()));
        list.addPolar(2, Math.PI);
        //Assertions.assertEquals(Complex.ofCis(pi), Complex.parse(Complex.ofCis(pi).toString()));
        list.addCis(Math.PI);

        return list;
    }

    /**
     * Test parse and toString are compatible.
     */
    @Test
    void testParseAndToString() {

        ComplexList list = createList();

        String listToString = list.toString();

        ComplexList parsedList = ComplexList.parse(listToString);

        Assertions.assertEquals(list, parsedList);

    }

    /**
     * Test list forEach conj.
     */
    @Test
    void testConj() {
        final ComplexList list = createList();

        final ComplexList originalCopy = ComplexList.parse(list.toString());

        //conjugate of conjugate - should get back original
        list.conj();

        Assertions.assertEquals(list, originalCopy);

        list.conj().multiply(Complex.ofCartesian(1, 1));
    }

    /**
     * Test list forEach conj.
     */
    @Test
    void testForEachConj() {
        final ComplexList list = createList();

        final ComplexList originalCopy = ComplexList.parse(list.toString());

        //conjugate of conjugate - should get back original
        list.forEach(ComplexFunctions::conj).forEach(ComplexFunctions::conj);

        Assertions.assertEquals(list, originalCopy);

        list.forEach(ComplexFunctions::conj).multiply(Complex.ofCartesian(1, 1));
    }
    /**
     * Test list forEach conj.
     */
    @Test
    void testArrayConj() {
        final ComplexList list = createList();

        final ComplexDoubleArray originalCopy = ComplexList.parse(list.toString());

        //conjugate of conjugate - should get back original
        list.apply(ComplexArrayFunctions::conj).apply(ComplexArrayFunctions::conj);

        Assertions.assertEquals(list, originalCopy);

    }
    /**
     * Test list multiply.
     */
    @Test
    void testMultiplyList() {
        final Complex x = Complex.ofCartesian(3.0, 4.0);
        final Complex y = Complex.ofCartesian(5.0, 6.0);
        final ComplexList list = new ComplexList();
        list.add(x);

        list.multiply(y);

        final ComplexDouble z = list.get(0);

        Assertions.assertEquals(-9.0, z.getReal());
        Assertions.assertEquals(38.0, z.getImaginary());
    }
    /**
     * Test list array exp
     */
    @Test
    void testArrayExp() {
        final ComplexList list = createList();

        final ComplexDoubleArray originalCopy = ComplexList.parse(list.toString());

        //exp


        Assertions.assertEquals(list, originalCopy);
    }
}
