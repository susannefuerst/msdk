package io.github.msdk.io.mzxml;

import java.io.File;
import java.util.List;
import java.util.Optional;

import javax.annotation.Nonnull;

import com.google.common.collect.ImmutableList;

import io.github.msdk.datamodel.chromatograms.Chromatogram;
import io.github.msdk.datamodel.files.FileType;
import io.github.msdk.datamodel.rawdata.MsFunction;
import io.github.msdk.datamodel.rawdata.MsScan;
import io.github.msdk.datamodel.rawdata.RawDataFile;

public class MzXMLRawDataFile implements RawDataFile {

  private static final @Nonnull FileType fileType = FileType.MZXML;

  private final @Nonnull File sourceFile;

  private final @Nonnull List<MsFunction> msFunctions;
  private final @Nonnull List<MsScan> msScans;
  private final @Nonnull List<Chromatogram> chromatograms;

  private @Nonnull String name;

  /**
   * <p>
   * Constructor for MzXMLRawDataFile.
   * </p>
   *
   * @param sourceFile a {@link java.io.File} object.
   * @param msFunctions a {@link java.util.List} object.
   * @param msScans a {@link java.util.List} object.
   * @param chromatograms a {@link java.util.List} object.
   */
  @SuppressWarnings("null")
  public MzXMLRawDataFile(@Nonnull File sourceFile, List<MsFunction> msFunctions,
      List<MsScan> msScans, List<Chromatogram> chromatograms) {
    this.sourceFile = sourceFile;
    this.name = sourceFile.getName();
    this.msFunctions = msFunctions;
    this.msScans = msScans;
    this.chromatograms = chromatograms;
  }

  /** {@inheritDoc} */
  @Override
  @Nonnull
  public String getName() {
    return name;
  }

  /** {@inheritDoc} */
  @Override
  @Nonnull
  public Optional<File> getOriginalFile() {
    return Optional.ofNullable(sourceFile);
  }

  /** {@inheritDoc} */
  @Override
  @Nonnull
  public FileType getRawDataFileType() {
    return fileType;
  }

  /** {@inheritDoc} */
  @SuppressWarnings("null")
  @Override
  @Nonnull
  public List<MsFunction> getMsFunctions() {
    return ImmutableList.copyOf(msFunctions);
  }

  /** {@inheritDoc} */
  @SuppressWarnings("null")
  @Override
  @Nonnull
  public List<MsScan> getScans() {
    return ImmutableList.copyOf(msScans);
  }

  /** {@inheritDoc} */
  @SuppressWarnings("null")
  @Override
  @Nonnull
  public List<Chromatogram> getChromatograms() {
    return ImmutableList.copyOf(chromatograms);
  }

  /** {@inheritDoc} */
  @Override
  public void dispose() {}

}
