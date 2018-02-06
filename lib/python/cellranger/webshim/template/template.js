var app = angular.module("webshim", ["ui.bootstrap"]);

function supportsWebGl() {
  try {
    var canvas = document.createElement('canvas');
    var ctx = canvas.getContext('webgl') || canvas.getContext('experimental-webgl');
    var exts = ctx.getSupportedExtensions();
  } catch (e) {
    return false;
  }

  return true;
};

function preProcess(data) {
  for (var i = 0; i < data.charts.length; i++) {
    var chart = data.charts[i];
    if (chart.layout) {
      chart.title = chart.layout.title;
      chart.layout.title = '';
    }
  }

  var supports_webgl = supportsWebGl();
  if (supports_webgl) {
    return data;
  }

  for (var i = 0; i < data.charts.length; i++) {
    var chart = data.charts[i];
    if (chart.data) {
      for (var j = 0; j < chart.data.length; j++) {
        var trace = chart.data[j];
        if (trace.type == "scattergl") {
          trace.type = "scatter";
        }
      }
    }
  }

  return data;
};

function show_description(ev) {
  var target = $(ev.target);
  target.closest(".has_desc").children(".summary_description").toggle();
};

app.filter('lineSplit', ['$sce', function($sce) {
    return function(input, splitChar) {
        return $sce.trustAsHtml(input.replace(new RegExp(splitChar, 'g'), '<br/>'));
    };
}]);

app.filter('percent1', function() {
    return function(input) {
        return (parseFloat(input) * 100).toFixed(1) + '%';
    }
});

app.controller('WebshimCtrl', function($scope) {
  $scope.data = preProcess(JSON.parse(LZString.decompressFromEncodedURIComponent(compressed_data)));
  $scope.charts = [];
  $scope.alerts = $scope.data.alarms;
  $scope.alerts_any_errors = _.find($scope.alerts, { level: 'error'} );
  $scope.alerts_show_all = false;
  $scope.ids = ["summary", "analysis"];
  $scope.tab = $scope.ids[0];
  $(".summary_description").hide();

  $scope.showAlerts = function() {
    $scope.alerts_show_all = !$scope.alerts_show_all;
    if ($scope.alerts_show_all) {
      $('#alerts').slideDown('fast');
    } else {
      $('#alerts').slideUp('fast');
    }
  };

  $scope.filterChart = function(chart) {
    if (!("filters" in chart)) {
      return true;
    }

    for (var title in $scope.data.filters) {
      if (!(title in chart.filters) ||
          chart.filters[title] == $scope.data.filters[title].selected) {
        return true;
      }
    }
    return false;
  };

  $scope.filterCharts = function() {
    $scope.charts = [];
    for (var i = 0; i < $scope.data.charts.length; i++) {
      var chart = $scope.data.charts[i];
      if ($scope.filterChart(chart)) {
        $scope.charts.push(chart);
      }
    }
  };

  $scope.selectFilterValue = function(filter, value) {
    filter.selected = value;
    $scope.filterCharts()
  };

  $scope.filterCharts();
});

app.directive("chartDiv", function() {
    return {
        scope: {
            chart: "="
        },
        link: function(scope, element, attrs) {
            scope.$watch(attrs.chartDiv, function(value) {
                if (scope.chart.table && scope.chart.name == "differential_expression") {
                    // init Sortable
                    $("#"+attrs.id).attr("data-sortable", "");
                    $("#"+attrs.id).addClass("sortable-theme-bootstrap");
                    scope.$evalAsync(function() {
                        Sortable.init();
                        $("#de_table_header-2").trigger("click");
                    });
                } else {
                    Plotly.newPlot(attrs.id, scope.chart.data, scope.chart.layout, scope.chart.config || {});
                }
            });
        }
    };
});
