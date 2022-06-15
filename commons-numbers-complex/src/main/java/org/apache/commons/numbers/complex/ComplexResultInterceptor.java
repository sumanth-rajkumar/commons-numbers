package org.apache.commons.numbers.complex;

import java.util.Stack;
import java.util.concurrent.atomic.AtomicLong;

public class ComplexResultInterceptor<R> implements ComplexResult<R> {

    private Stack<ComplexResult<?>> afterResultProviderStack = new Stack<>();
    private Stack<ComplexFunction<?>> afterFunctionStack= new Stack<>();
    private Stack<Object> afterResult = new Stack<>();

    public static final ThreadLocal<ComplexResultInterceptor> TLOCAL_ResultInterceptor = ThreadLocal.withInitial(() -> new ComplexResultInterceptor());

    private static final AtomicLong counter = new AtomicLong();

    private ComplexResultInterceptor() {
        System.out.println("Allocating ComplexResultInterceptor # " + counter.incrementAndGet());
    }

    public <U> void pushResultProvider(ComplexFunction<U> func, ComplexResult<U> provider) {
        afterFunctionStack.push(func);
        afterResultProviderStack.push(provider);
    }

    @Override
    public R apply(double r, double i) {
        ComplexFunction after = afterFunctionStack.pop();
        ComplexResult resultProvider = afterResultProviderStack.pop();
        afterResult.push(after.apply(r, i, resultProvider));
        return null;
    }

    public <U> U popResult() {
        return (U) afterResult.pop();
    }
}
