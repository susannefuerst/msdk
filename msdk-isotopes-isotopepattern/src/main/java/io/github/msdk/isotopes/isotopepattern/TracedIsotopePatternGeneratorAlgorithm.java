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
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map.Entry;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.annotation.Nonnull;

import org.openscience.cdk.config.IsotopeFactory;
import org.openscience.cdk.config.Isotopes;
import org.openscience.cdk.interfaces.IIsotope;

import com.google.common.base.Strings;

import io.github.msdk.datamodel.MsSpectrumType;
import io.github.msdk.isotopes.isotopepattern.impl.TracedIsotopePattern;
import io.github.msdk.util.MsSpectrumUtil;

/**
 * Class to simulate isotope patterns for tracing experiments. Assume you labeled with
 * (15N-amino)(13C)5-glutamine and you want to simulate the isotope pattern of the glutamate pool
 * where 60% of the pool consists of (15N)(13C)5-glutamate 0% of (15N)(12C)5-glutamate and 0% of
 * (14N)(13C)5-glutamate. Then your main parameters for the #simulateTracedPattern method are:
 * 
 * formula: C5H9NO4 (or C5H9NO4+) <br>
 * tracer1: 13C <br>
 * tracer2: 15N <br>
 * tracer1Inc: 0.0 <br>
 * tracer2Inc: 0.0 <br>
 * tracerAllInc: 0.6 <br>
 * 
 * 
 * If you labeled with (13C)6-glucose and (15N-amino)(12C)5-glutamine you may expect an independent
 * labeling resulting in (15N)(12C)5-glutamate, (15N)(13C)2(12C)3-glutamate and
 * (14N)(13C)2(12C)3-glutamate.To simulate a glutamate pool where 10% consist of
 * (15N)(12C)5-glutamate, 20% consist of (15N)(13C)2(12C)3-glutamate and 50% consist of
 * (14N)(13C)2(12C)3-glutamate then your main parameters for the #simulateTracedPattern method are:
 * 
 * formula: C5H9NO4 (or C5H9NO4+) <br>
 * tracer1: 13C <br>
 * tracer2: 15N <br>
 * tracer1Inc: 0.5 <br>
 * tracer2Inc: 0.1 <br>
 * tracerAllInc: 0.2 <br>
 * 
 * The input incorporation rates represent the real experimental incorporation (this is what you get
 * out of your spectrum after you corrected for natural abundances of isotopes). The simulated
 * spectrum itself is not corrected for natural abundances.
 * 
 */
public class TracedIsotopePatternGeneratorAlgorithm {


  private static final Pattern formulaPattern =
      Pattern.compile("^[\\[\\(]?(([A-Z][a-z]?[0-9]*)+)[\\]\\)]?(([0-9]*)([-+]))?$");
  private static final Pattern unchargedFormulaPattern =
      Pattern.compile("^(([A-Z][a-z]?[0-9]*)+)$");
  private static final Pattern tracerPattern = Pattern.compile("^([0-9]+)([A-Z][a-z]?)$");
  private static final Pattern multiElementPattern = Pattern.compile("([A-Z][a-z]?)([0-9]*)");
  private static final Pattern compositionPattern =
      Pattern.compile("\\[([0-9]+)\\]([A-Z][a-z]?)([0-9]*)");

  public static @Nonnull TracedIsotopePattern simulateTracedPattern(@Nonnull String chemicalFormula,
      @Nonnull String capacityFormula, String tracer1, String tracer2, @Nonnull double tracer1Inc,
      @Nonnull double tracer2Inc, @Nonnull double tracerAllInc, @Nonnull Double minAbundance,
      @Nonnull Float intensityScale, @Nonnull Double mzTolerance) throws IOException {
    parameterCheck(chemicalFormula, capacityFormula, tracer1, tracer2, tracer1Inc, tracer2Inc,
        tracerAllInc);
    Matcher m = formulaPattern.matcher(chemicalFormula);
    m.matches();
    String unchargedFormula = m.group(1);
    String chargeCount = m.group(4) == null ? "" : m.group(4);
    String chargeSign = m.group(5) == null ? "" : m.group(5);
    String chargeString = chargeCount + chargeSign;
    String tracer1ReducedFormula =
        "[" + reduceFormula(unchargedFormula, capacityFormula, tracer1, null) + "]" + chargeString;
    String tracer2ReducedFormula =
        "[" + reduceFormula(unchargedFormula, capacityFormula, null, tracer2) + "]" + chargeString;
    String tracerBothReducedFormula = "["
        + reduceFormula(unchargedFormula, capacityFormula, tracer1, tracer2) + "]" + chargeString;
    double totalInc = tracer1Inc + tracer2Inc + tracerAllInc;
    TracedIsotopePattern naturalPattern =
        (TracedIsotopePattern) IsotopePatternGeneratorAlgorithm.generateIsotopes(chemicalFormula,
            minAbundance, (float) (1.0f - totalInc), mzTolerance, true, false);
    TracedIsotopePattern tracer1ReducedPattern;
    if (tracer1ReducedFormula.contains("[]")) {
      tracer1ReducedPattern = new TracedIsotopePattern(new double[1],
          new float[] {(float) tracer1Inc}, 1, MsSpectrumType.CENTROIDED);
    } else {
      tracer1ReducedPattern =
          (TracedIsotopePattern) IsotopePatternGeneratorAlgorithm.generateIsotopes(
              tracer1ReducedFormula, minAbundance, (float) tracer1Inc, mzTolerance, true, false);
    }
    TracedIsotopePattern tracer2ReducedPattern;
    if (tracer2ReducedFormula.contains("[]")) {
      tracer2ReducedPattern = new TracedIsotopePattern(new double[1],
          new float[] {(float) tracer2Inc}, 1, MsSpectrumType.CENTROIDED);
    } else {
      tracer2ReducedPattern =
          (TracedIsotopePattern) IsotopePatternGeneratorAlgorithm.generateIsotopes(
              tracer2ReducedFormula, minAbundance, (float) tracer2Inc, mzTolerance, true, false);
    }
    TracedIsotopePattern tracerBothReducedPattern;
    if (tracerBothReducedFormula.contains("[]")) {
      tracerBothReducedPattern = new TracedIsotopePattern(new double[1],
          new float[] {(float) tracerAllInc}, 1, MsSpectrumType.CENTROIDED);
    } else {
      tracerBothReducedPattern = (TracedIsotopePattern) IsotopePatternGeneratorAlgorithm
          .generateIsotopes(tracerBothReducedFormula, minAbundance, (float) tracerAllInc,
              mzTolerance, true, false);
    }
    tracer1ReducedPattern = addTracerMass(tracer1ReducedPattern, capacityFormula, tracer1, null);
    tracer1ReducedPattern =
        addTracerComposition(tracer1ReducedPattern, capacityFormula, tracer1, null);
    tracer2ReducedPattern = addTracerMass(tracer2ReducedPattern, capacityFormula, null, tracer2);
    tracer2ReducedPattern =
        addTracerComposition(tracer2ReducedPattern, capacityFormula, null, tracer2);
    tracerBothReducedPattern =
        addTracerMass(tracerBothReducedPattern, capacityFormula, tracer1, tracer2);
    tracerBothReducedPattern =
        addTracerComposition(tracerBothReducedPattern, capacityFormula, tracer1, tracer2);
    TracedIsotopePattern mergedPattern = merge(naturalPattern, tracer1ReducedPattern);
    mergedPattern = merge(mergedPattern, tracer2ReducedPattern);
    mergedPattern = merge(mergedPattern, tracerBothReducedPattern);
    mergedPattern = normalize(mergedPattern, intensityScale);
    setHeavyIsotopes(mergedPattern);
    return mergedPattern;
  }

  private static void parameterCheck(@Nonnull String chemicalFormula,
      @Nonnull String capacityFormula, String tracer1, String tracer2, double tracer1Inc,
      double tracer2Inc, double tracerAllInc) {
    formulaCheck(chemicalFormula, formulaPattern);
    formulaCheck(capacityFormula, unchargedFormulaPattern);
    formulaCheck(tracer1, tracerPattern);
    formulaCheck(tracer2, tracerPattern);
    ratesCheck(tracer1Inc, tracer2Inc, tracerAllInc);
  }

  /**
   * Checks if sum of rates is in [0,1]
   * 
   * @param rates
   */
  private static void ratesCheck(double... rates) {
    double sum = 0.0;
    for (double rate : rates) {
      sum = sum + rate;
    }
    if (sum < 0 || sum > 1) {
      throw new IllegalArgumentException("(Sum of) rate(s) must be in [0,1].");
    }
  }

  /**
   * Checks if formula matches the formulaPattern.
   * 
   * @param formula
   * @param formulaPattern
   */
  private static void formulaCheck(String formula, Pattern formulaPattern) {
    if (formula == null) {
      return;
    }
    Matcher m = formulaPattern.matcher(formula);
    if (!m.matches())
      throw new IllegalArgumentException("Invalid formula: " + formula);
  }

  /**
   * 
   * @param formula expected format: C7H14NOSi
   * @return A {@link LinkedHashMap} that assigns to each element symbol its quantity given by the
   *         formula. e.g. [C:7, H:14, N:1, O:1, Si:1]
   */
  private static LinkedHashMap<String, Integer> formulaMap(String formula) {
    LinkedHashMap<String, Integer> elementFormula = new LinkedHashMap<String, Integer>();
    Matcher formulaMatcher = multiElementPattern.matcher(formula);
    ArrayList<String> elementTokens = new ArrayList<String>();
    while (formulaMatcher.find()) {
      elementTokens.add(formulaMatcher.group());
    }
    for (String elementToken : elementTokens) {
      Matcher elementMatcher = multiElementPattern.matcher(elementToken);
      if (elementMatcher.matches()) {
        String element = elementMatcher.group(1);
        Integer quantity = elementMatcher.group(2).equals("") ? Integer.valueOf(1)
            : Integer.valueOf(elementMatcher.group(2));
        elementFormula.put(element, quantity);
      }
    }
    return elementFormula;
  }

  /**
   * 
   * @param formula expected format [12]C7[1]H14[14]N[16]O[28]Si.
   * @return A {@link LinkedHashMap} that assigns to each isotope its quantity given by the formula.
   *         e.g. [12C:7, 1H:14, 14N:1, 16O:1, 28Si:1]
   */
  private static LinkedHashMap<String, Integer> isotopeMap(String formula) {
    LinkedHashMap<String, Integer> isotopeFormula = new LinkedHashMap<String, Integer>();
    Matcher formulaMatcher = compositionPattern.matcher(formula);
    ArrayList<String> isotopeTokens = new ArrayList<String>();
    while (formulaMatcher.find()) {
      isotopeTokens.add(formulaMatcher.group());
    }
    for (String isotopeToken : isotopeTokens) {
      Matcher isotopeMatcher = compositionPattern.matcher(isotopeToken);
      if (isotopeMatcher.matches()) {
        String isotope = isotopeMatcher.group(1) + isotopeMatcher.group(2);
        Integer quantity = isotopeMatcher.group(3).equals("") ? Integer.valueOf(1)
            : Integer.valueOf(isotopeMatcher.group(3));
        isotopeFormula.put(isotope, quantity);
      }
    }
    return isotopeFormula;
  }

  /**
   * 
   * @param formula expected format: C7H14NOSi.
   * @param capacity expected format: C4N.
   * @param tracer1 expected format: 13C
   * @param tracer2 expected format: 15N
   * @return The formula reduced by the tracer, according to their quantity given by the capacity
   *         formula. e.g. the input (C7H14NOSi, C4N, 13C, 15N) returns C3H14OSi; the input
   *         (C7H14NOSi, C4N, null, 15N) returns C7H14OSi
   */
  private static String reduceFormula(String formula, String capacity, String tracer1,
      String tracer2) {
    LinkedHashMap<String, Integer> elementFormula = formulaMap(formula);
    LinkedHashMap<String, Integer> capacityFormula = formulaMap(capacity);
    if (!Strings.isNullOrEmpty(tracer1)) {
      String tracer1Element = tracerElement(tracer1);
      if (elementFormula.get(tracer1Element) != null) {
        Integer quantity =
            capacityFormula.get(tracer1Element) != null ? capacityFormula.get(tracer1Element) : 0;
        elementFormula.put(tracer1Element, elementFormula.get(tracer1Element) - quantity);
      }
    }
    if (!Strings.isNullOrEmpty(tracer2)) {
      String tracer2Element = tracerElement(tracer2);
      if (elementFormula.get(tracer2Element) != null) {
        Integer quantity =
            capacityFormula.get(tracer2Element) != null ? capacityFormula.get(tracer2Element) : 0;
        elementFormula.put(tracer2Element, elementFormula.get(tracer2Element) - quantity);
      }
    }
    return toString(elementFormula);
  }

  /**
   * 
   * @param formula e.g. [C:7, H:14, N:1, O:1, Si:1]
   * @return Simple concatenation of all the key value pairs in formula. A key value pair like C:1
   *         will be reduced to only C (not C1). e.g C7H14NOSi.
   */
  private static String toString(LinkedHashMap<String, Integer> formula) {
    StringBuilder builder = new StringBuilder();
    for (Entry<String, Integer> entry : formula.entrySet()) {
      if (entry.getValue() != 0) {
        String number = entry.getValue() > 1 ? String.valueOf(entry.getValue()) : "";
        builder.append(entry.getKey() + number);
      }
    }
    return builder.toString();
  }

  /**
   * 
   * @param isotopeFormula e.g. [12C:7, 1H:14, 14N:1, 16O:1, 28Si:1]
   * @return A formatted concatenation of the key value pairs in isotopeFormula. e.g.
   *         [12]C7[1]H14[14]N[16]O[28]Si.
   */
  private static String toCompositionString(LinkedHashMap<String, Integer> isotopeFormula) {
    StringBuilder builder = new StringBuilder();
    for (Entry<String, Integer> entry : isotopeFormula.entrySet()) {
      if (entry.getValue() != null && entry.getValue() != 0) {
        String number = entry.getValue() > 1 ? String.valueOf(entry.getValue()) : "";
        builder.append(
            "[" + tracerMassNumber(entry.getKey()) + "]" + tracerElement(entry.getKey()) + number);
      }
    }
    return builder.toString();
  }

  /**
   * 
   * @param pattern
   * @param capacityFormula e.g. C4N
   * @param tracer1 e.g. 13C
   * @param tracer2 e.g. 15N
   * @return Adds the masses of tracer1 and tracer2 according to their quantity in capacityFormula
   *         to each mzValue of the pattern.
   * @throws IOException
   */
  private static TracedIsotopePattern addTracerMass(TracedIsotopePattern pattern,
      String capacityFormula, String tracer1, String tracer2) throws IOException {
    Double massToAdd = 0.0;
    if (tracer1 != null) {
      massToAdd = massToAdd + totalTracerMass(tracer1, capacityFormula);
    }
    if (tracer2 != null) {
      massToAdd = massToAdd + totalTracerMass(tracer2, capacityFormula);
    }
    double[] mzValues = pattern.getMzValues();
    for (int i = 0; i < mzValues.length; i++) {
      mzValues[i] = mzValues[i] + massToAdd;
    }
    String[] isotopeComposition = ((TracedIsotopePattern) pattern).getIsotopeComposition();
    return new TracedIsotopePattern(mzValues, pattern.getIntensityValues(), mzValues.length,
        pattern.getSpectrumType(), isotopeComposition);
  }

  /**
   * 
   * @param tracer e.g. 13C
   * @param capacityFormula e.g. C4N
   * @return the mass of the tracer multiplied by its quantity according to the capacityFormula. The
   *         input (13C, C4N) would return 4*massOf(13C).
   * @throws IOException
   */
  private static Double totalTracerMass(String tracer, String capacityFormula) throws IOException {
    IsotopeFactory isotopeFactory = Isotopes.getInstance();
    Double tracerMass =
        isotopeFactory.getIsotope(tracerElement(tracer), tracerMassNumber(tracer)).getExactMass();
    Integer factor = formulaMap(capacityFormula).get(tracerElement(tracer));
    factor = factor == null ? 0 : factor;
    return factor * tracerMass;
  }

  /**
   * 
   * @param tracer e.g 13C
   * @return The element symbol of the tracer.
   */
  private static String tracerElement(String tracer) {
    Matcher tracerMatcher = tracerPattern.matcher(tracer);
    tracerMatcher.matches();
    String tracerElement = tracerMatcher.group(2);
    return tracerElement != null ? tracerElement : "";
  }

  /**
   * 
   * @param tracer e.g. 13C
   * @return The mass number of the tracer.
   */
  private static int tracerMassNumber(String tracer) {
    Matcher tracerMatcher = tracerPattern.matcher(tracer);
    tracerMatcher.matches();
    String tracerMassNumber = tracerMatcher.group(1);
    return Integer.parseInt(tracerMassNumber);
  }

  /**
   * 
   * @param pattern1
   * @param pattern2
   * @return A {@link TracedIsotopePattern} that includes all the mzValues, intensities and
   *         isotopesCompositions of pattern1 and pattern2.
   */
  private static TracedIsotopePattern merge(TracedIsotopePattern pattern1,
      TracedIsotopePattern pattern2) {
    double[] mzValues1 = pattern1.getMzValues();
    float[] intensityValues1 = pattern1.getIntensityValues();
    double[] mzValues2 = pattern2.getMzValues();
    float[] intensityValues2 = pattern2.getIntensityValues();
    LinkedHashMap<Double, Float> massIntensityMap = new LinkedHashMap<Double, Float>();
    String[] compositions1 = pattern1.getIsotopeComposition();
    String[] compositions2 = pattern2.getIsotopeComposition();
    LinkedHashMap<Double, String> massCompositionMap = new LinkedHashMap<Double, String>();
    for (int i = 0; i < mzValues1.length; i++) {
      if (massIntensityMap.get(mzValues1[i]) != null) {
        massIntensityMap.put(mzValues1[i],
            massIntensityMap.get(mzValues1[i]) + intensityValues1[i]);
      } else {
        massIntensityMap.put(mzValues1[i], intensityValues1[i]);
      }
      massCompositionMap.put(mzValues1[i], compositions1[i]);
    }
    for (int i = 0; i < mzValues2.length; i++) {
      if (massIntensityMap.get(mzValues2[i]) != null) {
        massIntensityMap.put(mzValues2[i],
            massIntensityMap.get(mzValues2[i]) + intensityValues2[i]);
      } else {
        massIntensityMap.put(mzValues2[i], intensityValues2[i]);
      }
      massCompositionMap.put(mzValues2[i], compositions2[i]);
    }
    List<Entry<Double, Float>> massIntensityList = new ArrayList<>(massIntensityMap.entrySet());
    massIntensityList.sort(Entry.comparingByKey());
    int size = massIntensityList.size();
    double[] newMzValues = new double[size];
    float[] newIntensityValues = new float[size];
    String[] newCompositions = new String[size];
    for (int i = 0; i < size; i++) {
      newMzValues[i] = massIntensityList.get(i).getKey();
      newIntensityValues[i] = massIntensityList.get(i).getValue();
      newCompositions[i] = massCompositionMap.get(newMzValues[i]);
    }
    return new TracedIsotopePattern(newMzValues, newIntensityValues, size,
        pattern1.getSpectrumType(), newCompositions);
  }

  /**
   * 
   * @param pattern
   * @param intensityScale
   * @return A new {@link TracedIsotopePattern} with intensities of pattern normalized and scaled by
   *         intensityScale.
   */
  private static TracedIsotopePattern normalize(TracedIsotopePattern pattern,
      Float intensityScale) {
    float[] intensityValues = pattern.getIntensityValues();
    MsSpectrumUtil.normalizeIntensity(intensityValues, intensityValues.length, intensityScale);
    double[] mzValues = pattern.getMzValues();
    String[] isotopeComposition = pattern.getIsotopeComposition();
    return new TracedIsotopePattern(mzValues, intensityValues, mzValues.length,
        pattern.getSpectrumType(), isotopeComposition);
  }

  /**
   * 
   * @param pattern
   * @param capacityFormula e.g. C4N
   * @param tracer1 e.g. 13C
   * @param tracer2 e.g. 15N
   * @return A {@link TracedIsotopePattern} with all the mzValues, intensities and extended
   *         isotopeCompositions of pattern. Extended means the input (somePattern, C4N, 13C, 15N)
   *         would extend a composition string [12]C3[1]H14[16]O[28]Si to
   *         [12]C3[13]C4[1]H14[16]O[28]Si[15]N
   */
  private static TracedIsotopePattern addTracerComposition(TracedIsotopePattern pattern,
      String capacityFormula, String tracer1, String tracer2) {
    LinkedHashMap<String, Integer> capacity = formulaMap(capacityFormula);
    Integer tracer1Count = 0;
    if (tracer1 != null) {
      tracer1Count =
          capacity.get(tracerElement(tracer1)) != null ? capacity.get(tracerElement(tracer1)) : 0;
    }
    Integer tracer2Count = 0;
    if (tracer2 != null) {
      tracer2Count =
          capacity.get(tracerElement(tracer2)) != null ? capacity.get(tracerElement(tracer2)) : 0;
    }
    String[] compositions = pattern.getIsotopeComposition();
    for (int i = 0; i < compositions.length; i++) {
      LinkedHashMap<String, Integer> compositionMap = isotopeMap(compositions[i]);
      if (tracer1 != null) {
        if (compositionMap.get(tracer1) != null) {
          compositionMap.put(tracer1, compositionMap.get(tracer1) + tracer1Count);
        } else {
          compositionMap.put(tracer1, tracer1Count);
        }
      }
      if (tracer2 != null) {
        if (compositionMap.get(tracer2) != null) {
          compositionMap.put(tracer2, compositionMap.get(tracer2) + tracer2Count);
        } else {
          compositionMap.put(tracer2, tracer2Count);
        }
      }
      compositions[i] = toCompositionString(compositionMap);
    }
    pattern.setIsotopeComposition(compositions);
    return pattern;
  }

  private static void setHeavyIsotopes(TracedIsotopePattern pattern) throws IOException {
    IsotopeFactory isotopeFactory = Isotopes.getInstance();
    String[] heavyIsotopes = new String[pattern.getIsotopeComposition().length];
    String[] compositions = pattern.getIsotopeComposition();
    for (int i = 0; i < compositions.length; i++) {
      LinkedHashMap<String, Integer> compositionMap = isotopeMap(compositions[i]);
      LinkedHashMap<String, Integer> heavyIsotopesMap = new LinkedHashMap<String, Integer>();
      for (Entry<String, Integer> entry : compositionMap.entrySet()) {
        Integer massNumber = tracerMassNumber(entry.getKey());
        String elementSymbol = tracerElement(entry.getKey());
        IIsotope[] isotopes = isotopeFactory.getIsotopes(elementSymbol);
        Integer minimalIsotopeMass = 1000;
        for (IIsotope isotope : isotopes) {
          if (isotope.getNaturalAbundance() > 0) {
            minimalIsotopeMass = Math.min(minimalIsotopeMass, isotope.getMassNumber());
          }
        }
        if (massNumber > minimalIsotopeMass) {
          heavyIsotopesMap.put(entry.getKey(), entry.getValue());
        }
      }
      heavyIsotopes[i] = toCompositionString(heavyIsotopesMap);
    }
    pattern.setHeavyIsotopes(heavyIsotopes);
  }

}
