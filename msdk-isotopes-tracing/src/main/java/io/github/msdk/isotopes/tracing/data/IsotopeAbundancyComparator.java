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

package io.github.msdk.isotopes.tracing.data;

import java.util.Comparator;

import io.github.msdk.isotopes.tracing.data.constants.Isotope;

/**
 * Compares two isotopes by their abundance.
 * 
 * @author Susanne Fürst, susannefuerst@freenet.de, susanne.fuerst@mdc-berlin.de
 *
 */
public class IsotopeAbundancyComparator implements Comparator<Isotope> {

  /**
   * @return 0 if both have the same abundance, -1 if abundance of firstIsotope is less than
   *         abundance of otherIsotope, 1, otherwise.
   */
  @Override
  public int compare(Isotope firstIsotope, Isotope secondIsotope) {
    return firstIsotope.getAbundance().compareTo(secondIsotope.getAbundance());
  }

}
