/*
 * (C) Copyright 2015-2017 by MSDK Development Team
 *
 * This software is dual-licensed under either
 *
 * (a) the terms of the GNU Lesser General Public License version 2.1 as published by the Free
 * Software Foundation
 *
 * or (per the licensee's choosing)
 *
 * (b) the terms of the Eclipse Public License v1.0 as published by the Eclipse Foundation.
 */


package io.github.msdk.isotopes.isotopepattern;

import java.io.IOException;

import org.junit.Assert;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import io.github.msdk.isotopes.isotopepattern.impl.TracedIsotopePattern;

public class TracedIsotopePatternGeneratorAlgorithmTest {

  /**
   * For debugging change org.slf4j.simpleLogger.defaultLogLevel to "debug" in
   * simplelogger.properties.
   */
  public static final Logger LOG =
      LoggerFactory.getLogger(TracedIsotopePatternGeneratorAlgorithmTest.class);

  @Test
  public void test1() throws IOException {
    /*
     * Simulate a 13C tracing with an (natural abundance corrected) incorporation of 50%.
     */
    String chemicalFormula = "C";
    String capacityFormula = "C";
    String tracer1 = "13C";
    String tracer2 = null;
    double tracer1Inc = 0.5;
    double tracer2Inc = 0;
    double tracerAllInc = 0;
    double minAbundance = 0.0001;
    float intensityScale = 1f;
    double mzTolerance = 0.0001;
    TracedIsotopePattern pattern = TracedIsotopePatternGeneratorAlgorithm.simulateTracedPattern(
        chemicalFormula, capacityFormula, tracer1, tracer2, tracer1Inc, tracer2Inc, tracerAllInc,
        minAbundance, intensityScale, mzTolerance);
    pattern = pattern.roundMzValues(6);
    pattern = pattern.roundIntensities(6);
    pattern = pattern.intensitiesToMID();
    LOG.debug(pattern.toDebugString());
    Assert.assertEquals(new Integer(2), pattern.getNumberOfDataPoints());
    Assert.assertEquals(tracer1Inc, pattern.getIntensityValues()[1],
        /* allow for natural abundance uncertainty */0.01);
    Assert.assertEquals(1 - tracer1Inc, pattern.getIntensityValues()[0],
        /* allow for natural abundance uncertainty */0.01);
    Assert.assertEquals(12.0000, pattern.getMzValues()[0], 0.0001);
    Assert.assertEquals("", pattern.getHeavyIsotopes()[0]);
    Assert.assertEquals(13.0034, pattern.getMzValues()[1], 0.0001);
    Assert.assertEquals("[13]C", pattern.getHeavyIsotopes()[1]);
  }

  @Test
  public void test2() throws IOException {
    /*
     * Simulate a 13C tracing with an (natural abundance corrected) incorporation of 50%. Simulation
     * should ignore unnecessary 15N tracer, because there is no N capacity.
     */
    String chemicalFormula = "C";
    String capacityFormula = "C";
    String tracer1 = "13C";
    String tracer2 = "15N";
    double tracer1Inc = 0.5;
    double tracer2Inc = 0.5;
    double tracerAllInc = 0;
    double minAbundance = 0.0001;
    float intensityScale = 1f;
    double mzTolerance = 0.0001;
    TracedIsotopePattern pattern = TracedIsotopePatternGeneratorAlgorithm.simulateTracedPattern(
        chemicalFormula, capacityFormula, tracer1, tracer2, tracer1Inc, tracer2Inc, tracerAllInc,
        minAbundance, intensityScale, mzTolerance);
    pattern = pattern.roundMzValues(6);
    pattern = pattern.roundIntensities(6);
    pattern = pattern.intensitiesToMID();
    LOG.debug(pattern.toDebugString());
    Assert.assertEquals(new Integer(2), pattern.getNumberOfDataPoints());
    Assert.assertEquals(tracer1Inc, pattern.getIntensityValues()[1],
        /* allow for natural abundance uncertainty */0.01);
    Assert.assertEquals(1 - tracer1Inc, pattern.getIntensityValues()[0],
        /* allow for natural abundance uncertainty */0.01);
  }

  @Test
  public void test3() throws IOException {
    /*
     * Simulate a 13C tracing with an incorporation of 100%.
     */
    String chemicalFormula = "C";
    String capacityFormula = "C";
    String tracer1 = "13C";
    String tracer2 = null;
    double tracer1Inc = 1;
    double tracer2Inc = 0;
    double tracerAllInc = 0;
    double minAbundance = 0.0001;
    float intensityScale = 1f;
    double mzTolerance = 0.0001;
    TracedIsotopePattern pattern = TracedIsotopePatternGeneratorAlgorithm.simulateTracedPattern(
        chemicalFormula, capacityFormula, tracer1, tracer2, tracer1Inc, tracer2Inc, tracerAllInc,
        minAbundance, intensityScale, mzTolerance);
    pattern = pattern.roundMzValues(6);
    pattern = pattern.roundIntensities(6);
    pattern = pattern.intensitiesToMID();
    LOG.debug(pattern.toDebugString());
    Assert.assertEquals(new Integer(2), pattern.getNumberOfDataPoints());
    Assert.assertEquals(0, pattern.getIntensityValues()[0], 0);
    Assert.assertEquals(1, pattern.getIntensityValues()[1], 0);
  }

  @Test
  public void test4() throws IOException {
    /*
     * Simulate a 13C tracing with an (natural abundance corrected) incorporation of 50%. Simulation
     * should ignore unnecessary 15N tracer, because there is no N capacity.
     */
    String chemicalFormula = "CN";
    String capacityFormula = "C";
    String tracer1 = "13C";
    String tracer2 = "15N";
    double tracer1Inc = 0.5;
    double tracer2Inc = 0.5;
    double tracerAllInc = 0;
    double minAbundance = 0.0001;
    float intensityScale = 1f;
    double mzTolerance = 0.0001;
    TracedIsotopePattern pattern = TracedIsotopePatternGeneratorAlgorithm.simulateTracedPattern(
        chemicalFormula, capacityFormula, tracer1, tracer2, tracer1Inc, tracer2Inc, tracerAllInc,
        minAbundance, intensityScale, mzTolerance);
    pattern = pattern.roundMzValues(6);
    pattern = pattern.roundIntensities(6);
    pattern = pattern.intensitiesToMID();
    LOG.debug(pattern.toDebugString());
    Assert.assertEquals(new Integer(4), pattern.getNumberOfDataPoints());
    Assert.assertEquals(tracer1Inc, pattern.getIntensityValues()[2],
        /* allow for natural abundance uncertainty */0.01);
    Assert.assertEquals(1 - tracer1Inc, pattern.getIntensityValues()[0],
        /* allow for natural abundance uncertainty */0.01);
    Assert.assertEquals(26.0031, pattern.getMzValues()[0], 0.0001);
    Assert.assertEquals("", pattern.getHeavyIsotopes()[0]);
    Assert.assertEquals(27.0001, pattern.getMzValues()[1], 0.0001);
    Assert.assertEquals("[15]N", pattern.getHeavyIsotopes()[1]);
    Assert.assertEquals(27.0064, pattern.getMzValues()[2], 0.0001);
    Assert.assertEquals("[13]C", pattern.getHeavyIsotopes()[2]);
    Assert.assertEquals(28.0035, pattern.getMzValues()[3], 0.0001);
    Assert.assertEquals("[15]N[13]C", pattern.getHeavyIsotopes()[3]);
  }

  @Test
  public void test5() throws IOException {
    /*
     * Simulate a 15N tracing with an (natural abundance corrected) incorporation of 50%. Simulation
     * should ignore unnecessary 13C tracer, because there is no C capacity.
     */
    String chemicalFormula = "CN";
    String capacityFormula = "N";
    String tracer1 = "13C";
    String tracer2 = "15N";
    double tracer1Inc = 0.5;
    double tracer2Inc = 0.5;
    double tracerAllInc = 0;
    double minAbundance = 0.0001;
    float intensityScale = 1f;
    double mzTolerance = 0.0001;
    TracedIsotopePattern pattern = TracedIsotopePatternGeneratorAlgorithm.simulateTracedPattern(
        chemicalFormula, capacityFormula, tracer1, tracer2, tracer1Inc, tracer2Inc, tracerAllInc,
        minAbundance, intensityScale, mzTolerance);
    pattern = pattern.roundMzValues(6);
    pattern = pattern.roundIntensities(6);
    pattern = pattern.intensitiesToMID();
    LOG.debug(pattern.toDebugString());
    Assert.assertEquals(new Integer(4), pattern.getNumberOfDataPoints());
    Assert.assertEquals(tracer1Inc, pattern.getIntensityValues()[1],
        /* allow for natural abundance uncertainty */0.01);
    Assert.assertEquals(1 - tracer1Inc, pattern.getIntensityValues()[0],
        /* allow for natural abundance uncertainty */0.01);
  }

  @Test
  public void test6() throws IOException {
    /*
     * Simulate a 13C,15N tracing with an (natural abundance corrected) incorporation of 20% 13C
     * only, 20% 15N only and 20% 13C15N.
     */
    String chemicalFormula = "CN";
    String capacityFormula = "CN";
    String tracer1 = "13C";
    String tracer2 = "15N";
    double tracer1Inc = 0.2;
    double tracer2Inc = 0.2;
    double tracerAllInc = 0.2;
    double minAbundance = 0.0001;
    float intensityScale = 1f;
    double mzTolerance = 0.0001;
    TracedIsotopePattern pattern = TracedIsotopePatternGeneratorAlgorithm.simulateTracedPattern(
        chemicalFormula, capacityFormula, tracer1, tracer2, tracer1Inc, tracer2Inc, tracerAllInc,
        minAbundance, intensityScale, mzTolerance);
    pattern = pattern.roundMzValues(6);
    pattern = pattern.roundIntensities(6);
    pattern = pattern.intensitiesToMID();
    LOG.debug(pattern.toDebugString());
    Assert.assertEquals(new Integer(4), pattern.getNumberOfDataPoints());
    Assert.assertEquals(1 - tracer1Inc - tracer2Inc - tracerAllInc, pattern.getIntensityValues()[0],
        /* allow for natural abundance uncertainty */0.01);
    Assert.assertEquals(tracer2Inc, pattern.getIntensityValues()[1],
        /* allow for natural abundance uncertainty */0.01);
    Assert.assertEquals(tracer1Inc, pattern.getIntensityValues()[2],
        /* allow for natural abundance uncertainty */0.01);
    Assert.assertEquals(tracerAllInc, pattern.getIntensityValues()[3],
        /* allow for natural abundance uncertainty */0.01);
  }

  @Test
  public void test7() throws IOException {
    /*
     * Simulate a natural pattern.
     */
    String chemicalFormula = "CN";
    String capacityFormula = "CN";
    String tracer1 = null;
    String tracer2 = null;
    double tracer1Inc = 0;
    double tracer2Inc = 0;
    double tracerAllInc = 0;
    double minAbundance = 0.00001;
    float intensityScale = 1f;
    double mzTolerance = 0.0001;
    TracedIsotopePattern pattern = TracedIsotopePatternGeneratorAlgorithm.simulateTracedPattern(
        chemicalFormula, capacityFormula, tracer1, tracer2, tracer1Inc, tracer2Inc, tracerAllInc,
        minAbundance, intensityScale, mzTolerance);
    pattern = pattern.roundMzValues(6);
    pattern = pattern.roundIntensities(6);
    pattern = pattern.intensitiesToMID();
    LOG.debug(pattern.toDebugString());
    Assert.assertEquals(new Integer(4), pattern.getNumberOfDataPoints());
    Assert.assertEquals(0.9857, pattern.getIntensityValues()[0], 0.0001);
    Assert.assertEquals(0.0036, pattern.getIntensityValues()[1], 0.0001);
    Assert.assertEquals(0.0107, pattern.getIntensityValues()[2], 0.0001);
    Assert.assertEquals(0.0000, pattern.getIntensityValues()[3], 0.0001);
  }
}
