//! Functions and traits useful for testing the martian stages

use anyhow::Result;
use martian::prelude::*;
use std::path::PathBuf;

/// Create a directory within the rover, run the given stage in
/// that directory and return the stage outputs
pub fn run_stage<S>(
    stage: S,
    args: <S as MartianStage>::StageInputs,
    rover: &MartianRover,
) -> Result<<S as MartianStage>::StageOutputs>
where
    S: MartianStage + Sync,
    <S as MartianStage>::StageInputs: Clone + Send + Sync,
    <S as MartianStage>::ChunkInputs: Clone + Send + Sync,
    <S as MartianStage>::ChunkOutputs: Send + Sync,
{
    let run_dir: PathBuf = rover.make_path(S::stage_name());
    std::fs::create_dir(&run_dir)?;
    let outs = stage.test_run(&run_dir, args)?;
    Ok(outs)
}

/*

/// A trait useful for setting up a test for a stage from a ground truth pipestance.
/// You need to specify:
/// - Which stage you are testing
/// - Where the test artifacts will be stored
/// - How to setup args and expected outs for the stage from the ground truth pipestance
/// - Given one set of outputs, how do you check them against the expected outputs
///
/// You will get these through the default implementations:
/// - Setting up the test appropriately in the test folder given a name and a pipestance `setup_test()`
/// - Run the stage withouth checking the outputs `dry_run()`
/// - Run the stage and check the outputs `run_test()`
pub trait StageTesting
where
    Self::Stage: MartianStage + Sync,
    <Self::Stage as MartianStage>::StageInputs: CopyInto + Clone + Send + Sync + Serialize,
    <Self::Stage as MartianStage>::StageOutputs: CopyInto + DeserializeOwned,
    <Self::Stage as MartianStage>::ChunkInputs: Clone + Send + Sync,
    <Self::Stage as MartianStage>::ChunkOutputs: Send + Sync,
{
    type Stage;
    /// The folder where all the artifacts related to this test are stored
    fn test_base() -> PathBuf;

    /// Given a "ground truth" pipestance create the stage inputs.
    fn setup_args_from_pipestance(
        pipestance: &Path,
    ) -> Result<<Self::Stage as MartianStage>::StageInputs>;

    /// Given a "ground truth" pipestance create the expected stage outputs.
    fn setup_outs_from_pipestance(
        pipestance: &Path,
    ) -> Result<<Self::Stage as MartianStage>::StageOutputs>;

    fn check_correctness(
        actual_outs: <Self::Stage as MartianStage>::StageOutputs,
        expected_outs: <Self::Stage as MartianStage>::StageOutputs,
    ) -> Result<()>;

    fn stage() -> Self::Stage;

    /// The folder within the test base where we store files related to the
    /// test `test_name`
    fn test_folder(test_name: &Path) -> PathBuf {
        Self::test_base().join(test_name)
    }

    fn stage_args_file(test_name: &Path) -> JsonFile {
        JsonFile::new(Self::test_folder(test_name), "args")
    }

    fn stage_outs_file(test_name: &Path) -> JsonFile {
        JsonFile::new(Self::test_folder(test_name), "outs")
    }

    fn setup_test(test_name: &Path, pipestance: &Path) -> Result<()> {
        let test_base = Self::test_base();
        assert!(
            test_base.exists(),
            "{} directory does not exist. You need to create it",
            test_base.display()
        );

        let test_folder = Self::test_folder(test_name);
        assert!(
            !test_folder.exists(),
            "Test folder '{}' already exists. Please delete it if you want it regenerated",
            test_folder.display()
        );

        assert!(
            pipestance.exists(),
            "Pipestance {} does not exist. You need to supply the absolute path with no tilde(~)",
            pipestance.display()
        );

        std::fs::create_dir(&test_folder)?;
        set_permissions(&test_folder)?;

        // Create the stage inputs
        let args = Self::setup_args_from_pipestance(pipestance)?;
        // Copy the files within args to the test folder
        let inputs_folder = test_folder.join("inputs");
        std::fs::create_dir(&inputs_folder)?;
        set_permissions(&inputs_folder)?;
        let args = args.copy_into(&inputs_folder)?;
        // Write the args json
        let args_file = Self::stage_args_file(test_name);
        args_file.write(&args)?;
        set_permissions(&args_file)?;

        // Create the stage outputs
        let outs = Self::setup_outs_from_pipestance(pipestance)?;
        // Copy the files within outs to the test folder
        let outputs_folder = test_folder.join("expected_outs");
        std::fs::create_dir(&outputs_folder)?;
        set_permissions(&outputs_folder)?;
        let outs = outs.copy_into(&outputs_folder)?;
        // Write the outs json
        let outs_file = Self::stage_outs_file(test_name);
        outs_file.write(&outs)?;
        set_permissions(&outs_file)?;

        Ok(())
    }

    fn dry_run(
        test_name: &Path,
        run_dir: &Path,
    ) -> Result<<Self::Stage as MartianStage>::StageOutputs> {
        Self::stage().test_run(run_dir, Self::stage_args_file(test_name).read()?)
    }

    fn failed_folder(test_name: &Path) -> PathBuf {
        Self::test_folder(test_name).join(format!("failed-{test_name:?}"))
    }

    /// 1. Run the test in a temporary directory
    /// 2. Copy the outputs to a folder named `failed-xxxx`
    /// 3. Remove the temporary directory
    /// 4. Check the outputs against expected outputs
    /// 5. Remove the failed directory if the checks pass
    fn run_test(test_name: &Path) -> Result<()> {
        let run_dir = tempfile::tempdir()?;
        let actual_outs = Self::dry_run(test_name, run_dir.path())?;
        let expected_outs = Self::stage_outs_file(test_name).read()?;

        // Copy the outputs to the failed folder
        let fail_path = Self::failed_folder(test_name);
        std::fs::create_dir(&fail_path)?;
        set_permissions(&fail_path)?;
        let actual_outs = actual_outs.copy_into(&fail_path)?;
        println!("Copied the outputs to {fail_path:?}");
        drop(run_dir);

        // Check correctness
        Self::check_correctness(actual_outs, expected_outs)?;
        println!("Correctness checks pass!");

        // Sometimes the tests fail while removing the folder for unknown reasons.
        // If the attempt to remove fails, wait for a while and retry (upto 3 times)
        let mut num_attempts = 0;
        loop {
            num_attempts += 1;
            let result = std::fs::remove_dir_all(&fail_path)
                .with_context(|| format!("While removing {fail_path:?}"));
            match result {
                Ok(_) => break,
                Err(e) => {
                    if num_attempts >= 3 {
                        return Err(e);
                    }
                    // Wait for 5 seconds before attempting to remove the directory again
                    std::thread::sleep(Duration::from_secs(5));
                    println!("Removing directory {fail_path:?} failed. Retrying..");
                }
            }
        }
        Ok(())
    }
}

/// Set group write permission on coped files (linux only)
fn set_permissions(path: &Path) -> Result<()> {
    assert!(
        path.exists(),
        "Path {} does not exist to set any permission",
        path.display()
    );
    #[cfg(target_os = "linux")]
    {
        use std::os::unix::fs::PermissionsExt;
        let permissions = std::fs::Permissions::from_mode(0o775);
        std::fs::set_permissions(path, permissions)?;
    }
    Ok(())
}

*/
