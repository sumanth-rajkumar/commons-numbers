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
import java.util.function.Function;

@FunctionalInterface
public  interface ComplexResult<R> {

    R apply(double r, double i);

    default <V> ComplexResult<V> andThen(Function<? super R, ? extends V> after) {
        Objects.requireNonNull(after);
        return (r, i) -> after.apply(apply(r, i));
    }

    default ComplexResult<R> compose(ComplexFunction before) {
        Objects.requireNonNull(before);
        return (r, i) -> before.apply(r, i, (x, y) -> apply(x, y));
    }
}

