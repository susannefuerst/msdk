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

package io.github.msdk.isotopes.isotopepattern.impl;

import java.util.LinkedHashMap;
import java.util.Map.Entry;

import javax.annotation.Nonnull;

import io.github.msdk.datamodel.MsSpectrumType;
import io.github.msdk.datamodel.SimpleMsSpectrum;

public class TracedIsotopePattern extends SimpleMsSpectrum {

  private String[] isotopeComposition;
  private String[] heavyIsotopes;

  public TracedIsotopePattern(@Nonnull double[] mzValues, @Nonnull float[] intensityValues,
      @Nonnull Integer size, @Nonnull MsSpectrumType spectrumType, String[] isotopeComposition) {
    super(mzValues, intensityValues, size, spectrumType);
    this.isotopeComposition = isotopeComposition;
  }

  public TracedIsotopePattern(@Nonnull double[] mzValues, @Nonnull float[] intensityValues,
      @Nonnull Integer size, @Nonnull MsSpectrumType spectrumType, String[] isotopeComposition,
      String[] heavyIsotopes) {
    super(mzValues, intensityValues, size, spectrumType);
    this.isotopeComposition = isotopeComposition;
    this.heavyIsotopes = heavyIsotopes;
  }

  public TracedIsotopePattern(@Nonnull double[] mzValues, @Nonnull float[] intensityValues,
      @Nonnull Integer size, @Nonnull MsSpectrumType spectrumType) {
    super(mzValues, intensityValues, size, spectrumType);
    this.isotopeComposition = new String[] {""};
    this.heavyIsotopes = new String[] {""};
  }

  /**
   * 
   * @return The compositions of the peak inducing isotopologues.
   */
  public String[] getIsotopeComposition() {
    return isotopeComposition;
  }

  public void setIsotopeComposition(String[] isotopeCompostion) {
    this.isotopeComposition = isotopeCompostion;
  }

  /**
   * 
   * @return The peak inducing heavy isotopes.
   */
  public String[] getHeavyIsotopes() {
    return heavyIsotopes;
  }

  public void setHeavyIsotopes(String[] heavyIsotopes) {
    this.heavyIsotopes = heavyIsotopes;
  }

  public TracedIsotopePattern roundMzValues(int precision) {
    LinkedHashMap<Double, Float> massIntensityMap = new LinkedHashMap<Double, Float>();
    LinkedHashMap<Double, String> massCompositionMap = new LinkedHashMap<Double, String>();
    LinkedHashMap<Double, String> massHeavyIsotopesMap = new LinkedHashMap<Double, String>();
    double[] mzValues = getMzValues();
    float[] intensityValues = getIntensityValues();
    for (int i = 0; i < mzValues.length; i++) {
      double newMass =
          (double) ((Math.round(mzValues[i] * Math.pow(10, precision))) / Math.pow(10, precision));
      if (massIntensityMap.get(newMass) != null) {
        massIntensityMap.put(newMass, massIntensityMap.get(newMass) + intensityValues[i]);
      } else {
        massIntensityMap.put(newMass, intensityValues[i]);
      }
      massCompositionMap.put(newMass, isotopeComposition[i]);
      massHeavyIsotopesMap.put(newMass, heavyIsotopes[i]);
    }
    int newSize = massIntensityMap.size();
    double[] newMzValues = new double[newSize];
    float[] newIntensityValues = new float[newSize];
    int entryCount = 0;
    for (Entry<Double, Float> entry : massIntensityMap.entrySet()) {
      newMzValues[entryCount] = entry.getKey();
      newIntensityValues[entryCount] = entry.getValue();
      entryCount++;
    }
    entryCount = 0;
    String[] newCompositions = new String[newSize];
    for (Entry<Double, String> entry : massCompositionMap.entrySet()) {
      newCompositions[entryCount] = entry.getValue();
      entryCount++;
    }
    entryCount = 0;
    String[] newHeavyIsotopes = new String[newSize];
    for (Entry<Double, String> entry : massHeavyIsotopesMap.entrySet()) {
      newHeavyIsotopes[entryCount] = entry.getValue();
      entryCount++;
    }
    return new TracedIsotopePattern(newMzValues, newIntensityValues, newSize,
        this.getSpectrumType(), newCompositions, newHeavyIsotopes);
  }

  public TracedIsotopePattern intensitiesToMID() {
    float[] intensityValues = getIntensityValues();
    float sumOfIntensities = 0;
    for (int i = 0; i < intensityValues.length; i++) {
      sumOfIntensities = sumOfIntensities + intensityValues[i];
    }
    for (int i = 0; i < intensityValues.length; i++) {
      intensityValues[i] = intensityValues[i] / sumOfIntensities;
    }
    return new TracedIsotopePattern(getMzValues(), intensityValues, intensityValues.length,
        this.getSpectrumType(), isotopeComposition, heavyIsotopes);
  }

  public TracedIsotopePattern roundIntensities(int precision) {
    float[] intensityValues = getIntensityValues();
    for (int i = 0; i < intensityValues.length; i++) {
      intensityValues[i] = (float) ((Math.round(intensityValues[i] * Math.pow(10, precision)))
          / Math.pow(10, precision));
    }
    return new TracedIsotopePattern(getMzValues(), intensityValues, intensityValues.length,
        this.getSpectrumType(), isotopeComposition, heavyIsotopes);
  }

  public String toDebugString() {
    StringBuilder builder = new StringBuilder();
    builder.append("\n");
    builder.append(String.format("%-9s|", "Mass"));
    builder.append(String.format("%-9s|", "Intensity"));
    builder.append(String.format("%-20s|", "Heavy isotopes"));
    builder.append(String.format("%-50s|", "Isotope formula"));
    builder.append("\n");
    for (int i = 0; i < isotopeComposition.length; i++) {
      builder.append(String.format("%-5.6f|", getMzValues()[i]));
      builder.append(String.format("%-1.7f|", getIntensityValues()[i]));
      builder.append(String.format("%-20s|", heavyIsotopes[i]));
      builder.append(String.format("%-50s|", isotopeComposition[i]));
      builder.append("\n");
    }
    return builder.toString();
  }
}
