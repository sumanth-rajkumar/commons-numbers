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
package org.apache.commons.numbers.primes;

import java.text.MessageFormat;
import java.util.HashSet;
import java.util.List;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

class PrimesTest {

    static final int[] PRIMES = {
        //primes here have been verified one by one using Dario Alejandro Alpern's tool.
        //see http://www.alpertron.com.ar/ECM.HTM
        2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 43, 47, 53, 71, 73, 79, 89, 97,
        107, 137, 151, 157, 271, 293, 331, 409, 607, 617, 683, 829,
        1049, 1103, 1229, 1657,
        2039, 2053, //around first boundary in miller-rabin
        2251, 2389, 2473, 2699, 3271, 3389, 3449, 5653, 6449, 6869, 9067, 9091,
        11251, 12433, 12959, 22961, 41047, 46337, 65413, 80803, 91577, 92693,
        118423, 656519, 795659,
        1373639, 1373677, //around second boundary in miller-rabin
        588977, 952381,
        1013041, 1205999, 2814001,
        22605091,
        25325981, 25326023, //around third boundary in miller-rabin
        100000007, 715827881,
        2147483647//Integer.MAX_VALUE
        };

    static final int[] NOT_PRIMES = {
        //composite chosen at random + particular values used in algorithms such as boundaries for millerRabin
        4, 6, 8, 9, 10, 12, 14, 15, 16, 18, 20, 21, 22, 24, 25,
        275,
        2037, 2041, 2045, 2046, 2047, 2048, 2049, 2051, 2055, //around first boundary in miller-rabin
        9095,
        463465,
        1373637, 1373641, 1373651, 1373652, 1373653, 1373654, 1373655, 1373673, 1373675, 1373679, //around second boundary in miller-rabin
        25325979, 25325983, 25325993, 25325997, 25325999, 25326001, 25326003, 25326007, 25326009, 25326011, 25326021, 25326025, //around third boundary in miller-rabin
        100000005,
        1073741341, 1073741823, 2147473649, 2147483641, 2147483643, 2147483645, 2147483646};

    static final int[] BELOW_2 = {Integer.MIN_VALUE, -1, 0, 1};

    static final HashSet<Integer> PRIMES_SET = new HashSet<>();
    static {
        for (int p : PRIMES) {
            PRIMES_SET.add(p);
        }
    }

    void assertPrimeFactorsException(int n, String expected) {
        try {
            Primes.primeFactors(n);
            Assertions.fail("Exception not thrown");
        } catch (IllegalArgumentException e) {
            Assertions.assertEquals(expected, e.getMessage());
        }
    }

    void assertNextPrimeException(int n, String expected) {
        try {
            Primes.nextPrime(n);
            Assertions.fail("Exception not thrown");
        } catch (IllegalArgumentException e) {
            Assertions.assertEquals(expected, e.getMessage());
        }
    }

    @Test
    void testNextPrime() {

        Assertions.assertEquals(2, Primes.nextPrime(0));
        Assertions.assertEquals(2, Primes.nextPrime(1));
        Assertions.assertEquals(2, Primes.nextPrime(2));
        Assertions.assertEquals(3, Primes.nextPrime(3));
        Assertions.assertEquals(5, Primes.nextPrime(4));
        Assertions.assertEquals(5, Primes.nextPrime(5));

        for (int i = 0; i < SmallPrimes.PRIMES.length - 1; i++) {
            for (int j = SmallPrimes.PRIMES[i] + 1; j <= SmallPrimes.PRIMES[i + 1]; j++) {
                Assertions.assertEquals(SmallPrimes.PRIMES[i + 1], Primes.nextPrime(j));
            }
        }

        Assertions.assertEquals(25325981, Primes.nextPrime(25325981));
        for (int i = 25325981 + 1; i <= 25326023; i++) {
            Assertions.assertEquals(25326023, Primes.nextPrime(i));
        }

        Assertions.assertEquals(Integer.MAX_VALUE, Primes.nextPrime(Integer.MAX_VALUE - 10));
        Assertions.assertEquals(Integer.MAX_VALUE, Primes.nextPrime(Integer.MAX_VALUE - 1));
        Assertions.assertEquals(Integer.MAX_VALUE, Primes.nextPrime(Integer.MAX_VALUE));

        assertNextPrimeException(Integer.MIN_VALUE, MessageFormat.format(Primes.NUMBER_TOO_SMALL, Integer.MIN_VALUE, 0));
        assertNextPrimeException(-1, MessageFormat.format(Primes.NUMBER_TOO_SMALL, -1, 0));
        assertNextPrimeException(-13, MessageFormat.format(Primes.NUMBER_TOO_SMALL, -13, 0));
    }

    @Test
    void testIsPrime() throws Exception {
        for (int i : BELOW_2) {
            Assertions.assertFalse(Primes.isPrime(i));
        }
        for (int i:NOT_PRIMES) {
            Assertions.assertFalse(Primes.isPrime(i));
        }
        for (int i:PRIMES) {
            Assertions.assertTrue(Primes.isPrime(i));
        }
    }

    static int sum(List<Integer> numbers) {
        int out = 0;
        for (int i:numbers) {
            out += i;
        }
        return out;
    }

    static int product(List<Integer> numbers) {
        int out = 1;
        for (int i : numbers) {
            out *= i;
        }
        return out;
    }

    static void checkPrimeFactors(List<Integer> factors) {
        for (int p : factors) {
            if (!PRIMES_SET.contains(p)) {
                Assertions.fail("Not found in primes list: " + p);
            }
        }
    }

    @Test
    void testPrimeFactors() throws Exception {
        for (int i : BELOW_2) {
            assertPrimeFactorsException(i, MessageFormat.format(Primes.NUMBER_TOO_SMALL, i, 2));
        }
        for (int i : NOT_PRIMES) {
            List<Integer> factors = Primes.primeFactors(i);
            checkPrimeFactors(factors);
            int prod = product(factors);
            Assertions.assertEquals(i, prod);
        }
        for (int i : PRIMES) {
            List<Integer> factors = Primes.primeFactors(i);
            Assertions.assertEquals(i, (int)factors.get(0));
            Assertions.assertEquals(1, factors.size());
        }
    }
}
