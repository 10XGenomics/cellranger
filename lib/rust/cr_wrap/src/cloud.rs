use crate::env::PkgEnv;
use crate::utils::AllArgs;
use anyhow::Result;
use serde_json;
use std::env;
use std::fs::read_to_string;
use std::path::Path;
use std::process::ExitCode;
use url::Url;

pub fn run_cloud(pkg_env: &PkgEnv, args: &AllArgs) -> Result<ExitCode> {
    // target the test-cloud environment if the right environment variables are set
    let test_server = env::var("TENX_CLOUD_URL").unwrap_or_default();
    let tenx_cloud_token_path = env::var("TENX_CLOUD_TOKEN_PATH").unwrap_or_default();
    let cloud_cmd = "bin/tenkit/txg";
    if !test_server.is_empty() && !tenx_cloud_token_path.is_empty() {
        let mut new_input = args.input.clone();
        // read CF access token and add token params as headers to the request
        if let Ok((client_id, client_secret)) =
            extract_cloud_headers(&test_server, Path::new(&tenx_cloud_token_path))
        {
            new_input.push(String::from("-H"));
            new_input.push(format!("cf-access-client-id: {}", client_id.clone()));
            new_input.push(String::from("-H"));
            new_input.push(format!(
                "cf-access-client-secret: {}",
                client_secret.clone()
            ));
        } else {
            eprintln!("Could not find token path information at {tenx_cloud_token_path}");
        }

        new_input.push(String::from("--debug-url"));
        new_input.push(test_server);

        let new_args = &AllArgs { input: new_input };
        return pkg_env.run_subcmd(cloud_cmd, new_args);
    }
    pkg_env.run_subcmd(cloud_cmd, args)
}

// Extract a cf-access-client-id and cf-access-client-secret from the 10xcloud_token.json
// file located at the specified path.
fn extract_cloud_headers(server: &str, tenx_cloud_token_path: &Path) -> Result<(String, String)> {
    // assume debug url has complete scheme and URL
    let server_url = Url::parse(server)?;
    let server_host = server_url.host_str().expect("No host in token url");
    let token_contents = read_to_string(tenx_cloud_token_path)?;
    let token_dict: serde_json::Value = serde_json::from_str(&token_contents)?;
    let server_dict = &token_dict[server_host];
    // could do serde deserialization for this object but meh
    let cf_dict = &server_dict["cloudflare_token"];
    let client_id = &cf_dict["cf-access-client-id"];
    let client_secret = &cf_dict["cf-access-client-secret"];

    if client_id.is_string() && client_secret.is_string() {
        // to_string returns the quotes, weird
        return Ok((
            String::from(client_id.as_str().unwrap()),
            String::from(client_secret.as_str().unwrap()),
        ));
    }
    Err(anyhow::anyhow!("Could not read CloudFlare token"))
}
