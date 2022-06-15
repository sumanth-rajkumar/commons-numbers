package org.apache.commons.numbers.complex;

@FunctionalInterface
public interface DoubleComplexBiFunction {
    <R> R apply(double r, double i, double f, ComplexResult<R> result);
}
