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

    /**
     * Test parse and toString are compatible.
     */
    @Test
    void testParseAndToString() {
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

        String listToString = list.toString();

        ComplexList parsedList = ComplexList.parse(listToString);

        Assertions.assertEquals(list, parsedList);
    }
}
