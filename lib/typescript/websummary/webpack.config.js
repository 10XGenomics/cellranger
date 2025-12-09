const path = require("path");
const MiniCssExtractPlugin = require("mini-css-extract-plugin");
const CssMinimizerPlugin = require("css-minimizer-webpack-plugin");
const HtmlWebpackPlugin = require("html-webpack-plugin");

// Assume production by default so that development is opt-in.
const isProduction = process.env.NODE_ENV !== "development";

const devPlugins = [new HtmlWebpackPlugin({ template: "public/index.html" })];

const productionPlugins = [
  new MiniCssExtractPlugin({ filename: "tenx-websummary-styles.min.css" }),
];

const plugins = isProduction
  ? productionPlugins
  : [...devPlugins, ...productionPlugins];

module.exports = {
  devServer: {
    historyApiFallback: true,
    hot: true,
    open: true,
    port: 3000,
    static: {
      // Serve files from the 'public' directory so that fetches on / read into
      // public/
      directory: path.join(__dirname, "public"),
    },
    // Watch for changes in data.json. Note this is gitignored, so you'll have
    // to bring your own.
    watchFiles: ["public/data.json"],
  },
  // Ensures we have source maps for debugging the output in a browser.
  devtool: isProduction ? "source-map" : "eval-source-map",
  entry: "./index.ts",
  mode: isProduction ? "production" : "development",
  module: {
    rules: [
      {
        exclude: /node_modules\/(?!websummary-test)/,
        loader: "babel-loader",
        test: /\.js$/,
      },
      // for plotly https://github.com/plotly/plotly-webpack#explanations
      { test: /\.js$/, use: ["ify-loader"] },
      // Mostly copied from the docs; removed style-loader because we dont want
      // CSS in the page and instead use ExtractTextPlugin to generate a file
      // https://getbootstrap.com/docs/4.0/getting-started/webpack/#importing-precompiled-sass
      {
        test: /\.(scss)$/,
        use: [
          { loader: MiniCssExtractPlugin.loader },
          // translates CSS into CommonJS modules
          { loader: "css-loader" },
          {
            loader: "postcss-loader", // Run post css actions
            options: {
              // post css plugins, can be exported to postcss.config.js
              postcssOptions: {
                plugins: ["autoprefixer", "postcss-nested"],
              },
            },
          },
          // compiles Sass to CSS
          { loader: "sass-loader" },
        ],
      },
      {
        exclude: /node_modules|\.test\.tsx?$/,
        test: /\.tsx?$/,
        use: ["ts-loader"],
      },
      {
        test: /\.css$/,
        use: ["style-loader", "css-loader"],
      },
      {
        test: /\.(jpe?g|png|ttf|eot|svg|woff(2)?|ico)$/,
        type: "asset/inline",
      },
    ],
  },
  optimization: {
    minimize: isProduction,
    minimizer: [
      // Note:
      // Trying to do JS minimization in webpack in production should be avoided
      // due to memory consumption. A good sign that we're doing this wrong is
      // webpack failing during Bazel runs.
      new CssMinimizerPlugin({
        minimizerOptions: {
          preset: [
            "default",
            {
              discardComments: { removeAll: true },
            },
            {
              autoprefixer: {
                browsers: [
                  "chrome 79",
                  "edge 94",
                  "firefox 78",
                  "opera 80",
                  "safari 13.1",
                ],
              },
            },
          ],
        },
      }),
    ],
  },
  output: {
    // Clean the output directory before building for production
    clean: isProduction,
    filename: isProduction
      ? "tenx-websummary-script.min.js"
      : "tenx-websummary-script.js",
    path: path.resolve(__dirname, "./dist"),
  },
  plugins,
  // When building in a bazel sandbox, everything is a symlink.  Resolving
  // those symlinks breaks everything.
  resolve: {
    extensions: [".tsx", ".ts", ".js"],
    modules: [path.resolve(__dirname, "node_modules")],
    symlinks: false,
  },
};
