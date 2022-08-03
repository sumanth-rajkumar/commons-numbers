package org.apache.commons.numbers.complex;

import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.ListIterator;

public class ComplexArrayList implements List<ComplexDouble> {

    protected double[] realAndImagParts;

    private int size;

    public ComplexArrayList(double[] realAndImagParts, int index, int length) {
        this.realAndImagParts = realAndImagParts;
        this.size = length - index;
    }

    @Override
    public int size() {
        return size;
    }

    @Override
    public boolean isEmpty() {
        return size == 0;
    }

    @Override
    public boolean contains(Object o) {
        throw new UnsupportedOperationException();
    }

    @Override
    public Iterator<ComplexDouble> iterator() {
        throw new UnsupportedOperationException();
    }

    @Override
    public Object[] toArray() {
        throw new UnsupportedOperationException();
    }

    @Override
    public <T> T[] toArray(T[] a) {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean add(ComplexDouble complexDouble) {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean remove(Object o) {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean containsAll(Collection<?> c) {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean addAll(Collection<? extends ComplexDouble> c) {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean addAll(int index, Collection<? extends ComplexDouble> c) {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean removeAll(Collection<?> c) {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean retainAll(Collection<?> c) {
        throw new UnsupportedOperationException();
    }

    @Override
    public void clear() {
        throw new UnsupportedOperationException();
    }

    @Override
    public ComplexDouble get(int index) {
        return realAndImagParts[index];
    }

    @Override
    public ComplexDouble set(int index, ComplexDouble element) {
        return null;
    }

    @Override
    public void add(int index, ComplexDouble element) {
        throw new UnsupportedOperationException();
    }

    @Override
    public ComplexDouble remove(int index) {
        throw new UnsupportedOperationException();
    }

    @Override
    public int indexOf(Object o) {
        throw new UnsupportedOperationException();
    }

    @Override
    public int lastIndexOf(Object o) {
        throw new UnsupportedOperationException();
    }

    @Override
    public ListIterator<ComplexDouble> listIterator() {
        throw new UnsupportedOperationException();
    }

    @Override
    public ListIterator<ComplexDouble> listIterator(int index) {
        throw new UnsupportedOperationException();
    }

    @Override
    public List<ComplexDouble> subList(int fromIndex, int toIndex) {
        throw new UnsupportedOperationException();
    }
}
