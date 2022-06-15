package org.apache.commons.numbers.complex;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import java.util.concurrent.CompletableFuture;

public class ComplexComposeTest {


    private static ComplexFunction<Complex> neg = ComplexFunctions::negate;
    private static ComplexFunction<Complex> multiplyImag = ComplexFunctions::multiplyIma;
    private static ComplexFunction<Complex> conj = ComplexFunctions::conj;
    private static ComplexFunction<Complex> multiplyImagConj = multiplyImag.thenApply(conj);
    private static ComplexFunction<Complex> conjMultiplyImag = conj.thenApply(multiplyImag);
    private static ComplexFunction<Complex> identity1 = multiplyImagConj.thenApply(multiplyImagConj);
    private static ComplexFunction<Complex> identity2 = conjMultiplyImag.thenApply(conjMultiplyImag);






    public static void main(String[] args) throws Exception {
        int concurrency = 10;
        CompletableFuture[] futures = new CompletableFuture[concurrency];
        for (int i = 0; i < concurrency; i++) {
            futures[i] = CompletableFuture.runAsync(ComplexComposeTest::run);
        }
        CompletableFuture.allOf(futures).get();
    }

    @Test
    static void run() {
        double real = 1;
        double imag = 2;

        Complex c = Complex.ofCartesian(real, imag);

        System.out.println("real= " + real + " imag= " + imag);

        Complex c1 = multiplyImagConj.thenApply(neg).apply(c, Complex::ofCartesian);
        Complex c2 = conjMultiplyImag.apply(c, Complex::ofCartesian);
          


      
                  Assertions.assertEquals(c1, c2);
                              System.out.println("c1= " + c1 + ": c2= " + c2); 
        c1 = identity1.apply(c, Complex::ofCartesian);
        c2 = identity2.apply(c, Complex::ofCartesian);
                                 System.out.println("c1= " + c1 + ": c2= " + c2);

                           Assertions.assertEquals(c, c1);
                           Assertions.assertEquals(c, c2);

       
        
        System.out.println("c1= " + c1 + ": c2= " + c2);


    }
}
