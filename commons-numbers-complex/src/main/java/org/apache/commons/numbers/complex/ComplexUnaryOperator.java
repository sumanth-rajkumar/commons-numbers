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

import java.util.Objects;

/**
 * Represents a complex operation that accepts a complex number of type ComplexDouble and produces a ComplexDouble result.
 * @param <R> Generic
 */
@FunctionalInterface
public interface ComplexUnaryOperator<R> {

    /**
     * Represents an operator that accepts a complex number and a complex constructor to produce and return the result.
     * @param in Complex number
     * @param out Constructor
     * @return ComplexDouble
     */
    default R apply(ComplexDouble in, ComplexConstructor<R> out) {
        return apply(in.getReal(), in.getImaginary(), out);
    }

    /**
     * Represents an operator that accepts real and imaginary parts of a complex number and a complex constructor to produce and return the result.
     * @param r real
     * @param i imaginary
     * @param out Constructor
     * @return ComplexDouble
     */
    R apply(double r, double i, ComplexConstructor<R> out);

    /**
     * Returns a composed unary operator that first applies this function to its input, and then applies the after UnaryOperator function to the result.
     * If evaluation of either function throws an exception, it is relayed to the caller of the composed function.
     * @param after the function to apply after this function is applied
     * @return a composed function that first applies this function and then applies the after function
     */
    default ComplexUnaryOperator<R> andThen(ComplexUnaryOperator<R> after) {
        Objects.requireNonNull(after);
        return (real, imaginary, out) -> apply(real, imaginary, out.compose(after));
    }
}
