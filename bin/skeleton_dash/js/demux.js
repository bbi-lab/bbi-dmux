/********** Define components **********/
function RTTab(props) {
  var bad_wells_exist = Object.entries(props.bad_wells_barcodes).length === 0 && props.bad_wells_barcodes.constructor === Object;
  var plates = plate_list.map(function (plate_name, index) {
    return React.createElement(RTTabPlate, { key: index, lane: props.lane, plate_name: plate_name, norm: props.norm });
  });
  return React.createElement(
    "div",
    { className: props.className, id: "navrt-lane" + props.lane, role: "tabpanel", "aria-labelledby": "navrt-lane" + props.lane + "-tab" },
    plates,
    bad_wells_exist ? React.createElement("span", null) : React.createElement(BadWellTable, { bad_wells: props.bad_wells, lane: props.lane, bad_wells_barcodes: props.bad_wells_barcodes })
  );
}

function RTTabPlate(props) {
  var norm = props.norm;
  return React.createElement(
    "span",
    null,
    React.createElement(
      "h4",
      null,
      "Plate ",
      props.plate_name
    ),
    React.createElement("img", { src: "img/L00" + props.lane + "_" + props.plate_name + ".rt_plate.png", width: "50%", className: "rounded mx-auto d-block", alt: "..." }),
    norm ? React.createElement(
      "span",
      null,
      React.createElement(
        "h5",
        null,
        "Sentinel Normalized"
      ),
      React.createElement("img", { src: "img/L00" + props.lane + "_" + props.plate_name + ".rt_plate_sent_norm.png", width: "50%", className: "rounded mx-auto d-block", alt: "..." }),
      React.createElement(
        "h5",
        null,
        "Barnyard Normalized"
      ),
      React.createElement("img", { src: "img/L00" + props.lane + "_" + props.plate_name + ".rt_plate_barn_norm.png", width: "50%", className: "rounded mx-auto d-block", alt: "..." })
    ) : React.createElement("span", null)
  );
}

function TitleRow(props) {
  return React.createElement(
    "th",
    { scope: "col" },
    props.samp
  );
}

function RegRow(props) {
  return React.createElement(
    "td",
    null,
    props.val
  );
}

function NamedRow(props) {
  return React.createElement(
    "tr",
    null,
    React.createElement(
      "th",
      { scope: "row" },
      props.name
    ),
    React.createElement(
      "td",
      null,
      props.val
    )
  );
}

function BadWellTable(props) {
  var Lane = "Lane " + props.lane;
  return React.createElement(
    "span",
    null,
    React.createElement(
      "h4",
      null,
      "Wells outside plate layout for " + Lane
    ),
    React.createElement(
      "table",
      { className: "table table-hover" },
      React.createElement(
        "thead",
        null,
        React.createElement(
          "tr",
          null,
          React.createElement("th", { scope: "col" }),
          React.createElement(
            "th",
            { scope: "col" },
            "ReadCount"
          )
        )
      ),
      React.createElement(
        "tbody",
        null,
        props.bad_wells_barcodes.map(function (item, index) {
          return React.createElement(NamedRow, { key: index, name: item, val: props.bad_wells[[Lane]][[item]] });
        })
      )
    )
  );
}

function RTBarcodes(props) {
  return React.createElement(
    "div",
    { className: "tab-pane fade", id: "rt", role: "tabpanel", "aria-labelledby": "rt-tab" },
    React.createElement(
      "div",
      { className: "d-flex justify-content-between flex-wrap flex-md-nowrap align-items-center pt-3 pb-2 mb-3 border-bottom" },
      React.createElement(
        "h1",
        { className: "h3", id: "rt" },
        "RT Barcodes"
      )
    ),
    React.createElement(
      "nav",
      null,
      React.createElement(
        "div",
        { className: "nav nav-tabs", id: "navrt-tab", role: "tablist" },
        React.createElement(
          "a",
          { className: "nav-item nav-link", id: "navrt-lane1-tab", "data-toggle": "tab", href: "#navrt-lane1", role: "tab", "aria-controls": "navrt-lane1", "aria-selected": "false" },
          "Lane 1"
        ),
        React.createElement(
          "a",
          { className: "nav-item nav-link", id: "navrt-lane2-tab", "data-toggle": "tab", href: "#navrt-lane2", role: "tab", "aria-controls": "navrt-lane2", "aria-selected": "false" },
          "Lane 2"
        ),
        React.createElement(
          "a",
          { className: "nav-item nav-link", id: "navrt-lane3-tab", "data-toggle": "tab", href: "#navrt-lane3", role: "tab", "aria-controls": "navrt-lane3", "aria-selected": "false" },
          "Lane 3"
        ),
        React.createElement(
          "a",
          { className: "nav-item nav-link", id: "navrt-lane4-tab", "data-toggle": "tab", href: "#navrt-lane4", role: "tab", "aria-controls": "navrt-lane4", "aria-selected": "false" },
          "Lane 4"
        )
      )
    ),
    React.createElement(
      "div",
      { className: "tab-content", id: "nav-rtContent" },
      props.rt_tabs
    )
  );
}

function PCRTab(props) {
  var plates = pcr_combo_list.map(function (plate_name, index) {
    return React.createElement(PCRTabPlate, { key: index, lane: props.lane, plate_name: plate_name });
  });
  return React.createElement(
    "div",
    { className: props.className, id: "navpcr-lane" + props.lane, role: "tabpanel", "aria-labelledby": "navpcr-lane" + props.lane + "-tab" },
    React.createElement(
      "p",
      null,
      "Wells circled in ",
      React.createElement(
        "strong",
        null,
        React.createElement(
          "span",
          { style: { color: "red" } },
          "red"
        )
      ),
      " have >2 orders of magnitude fewer reads than the median well in the plate."
    ),
    plates
  );
}

function PCRTabPlate(props) {
  return React.createElement(
    "span",
    null,
    React.createElement(
      "h4",
      null,
      "PCR Combo ",
      props.plate_name
    ),
    React.createElement("img", { src: "img/L00" + props.lane + "_" + props.plate_name + ".pcr_plate.png", width: "50%", className: "rounded mx-auto d-block", alt: "..." })
  );
}

function PCRBarcodes(props) {
  return React.createElement(
    "div",
    { className: "tab-pane fade", id: "pcr", role: "tabpanel", "aria-labelledby": "pcr-tab" },
    React.createElement(
      "div",
      { className: "d-flex justify-content-between flex-wrap flex-md-nowrap align-items-center pt-3 pb-2 mb-3 border-bottom" },
      React.createElement(
        "h1",
        { className: "h3", id: "pcr" },
        "PCR Barcodes"
      )
    ),
    React.createElement(
      "nav",
      null,
      React.createElement(
        "div",
        { className: "nav nav-tabs", id: "navpcr-tab", role: "tablist" },
        React.createElement(
          "a",
          { className: "nav-item nav-link", id: "navpcr-lane1-tab", "data-toggle": "tab", href: "#navpcr-lane1", role: "tab", "aria-controls": "navpcr-lane1", "aria-selected": "false" },
          "Lane 1"
        ),
        React.createElement(
          "a",
          { className: "nav-item nav-link", id: "navpcr-lane2-tab", "data-toggle": "tab", href: "#navpcr-lane2", role: "tab", "aria-controls": "navpcr-lane2", "aria-selected": "false" },
          "Lane 2"
        ),
        React.createElement(
          "a",
          { className: "nav-item nav-link", id: "navpcr-lane3-tab", "data-toggle": "tab", href: "#navpcr-lane3", role: "tab", "aria-controls": "navpcr-lane3", "aria-selected": "false" },
          "Lane 3"
        ),
        React.createElement(
          "a",
          { className: "nav-item nav-link", id: "navpcr-lane4-tab", "data-toggle": "tab", href: "#navpcr-lane4", role: "tab", "aria-controls": "navpcr-lane4", "aria-selected": "false" },
          "Lane 4"
        )
      )
    ),
    React.createElement(
      "div",
      { className: "tab-content", id: "nav-pcrContent" },
      props.pcr_tabs
    )
  );
}

function LigTab(props) {
  var plates = lig_combo_list.map(function (plate_name, index) {
    return React.createElement(LigTabPlate, { key: index, lane: props.lane, plate_name: plate_name });
  });
  return React.createElement(
    "div",
    { className: props.className, id: "navlig-lane" + props.lane, role: "tabpanel", "aria-labelledby": "navlig-lane" + props.lane + "-tab" },
    plates
  );
}

function LigTabPlate(props) {
  return React.createElement(
    "span",
    null,
    React.createElement(
      "h4",
      null,
      "Ligation Plate ",
      props.plate_name
    ),
    React.createElement("img", { src: "img/L00" + props.lane + "_" + props.plate_name + ".lig_plate.png", width: "50%", className: "rounded mx-auto d-block", alt: "..." })
  );
}

function LigBarcodes(props) {
  return React.createElement(
    "div",
    { className: "tab-pane fade", id: "lig", role: "tabpanel", "aria-labelledby": "lig-tab" },
    React.createElement(
      "div",
      { className: "d-flex justify-content-between flex-wrap flex-md-nowrap align-items-center pt-3 pb-2 mb-3 border-bottom" },
      React.createElement(
        "h1",
        { className: "h3", id: "lig-name" },
        "Ligation Barcodes"
      )
    ),
    React.createElement(
      "nav",
      null,
      React.createElement(
        "div",
        { className: "nav nav-tabs", id: "navlig-tab", role: "tablist" },
        React.createElement(
          "a",
          { className: "nav-item nav-link", id: "navlig-lane1-tab", "data-toggle": "tab", href: "#navlig-lane1", role: "tab", "aria-controls": "navlig-lane1", "aria-selected": "false" },
          "Lane 1"
        ),
        React.createElement(
          "a",
          { className: "nav-item nav-link", id: "navlig-lane2-tab", "data-toggle": "tab", href: "#navlig-lane2", role: "tab", "aria-controls": "navlig-lane2", "aria-selected": "false" },
          "Lane 2"
        ),
        React.createElement(
          "a",
          { className: "nav-item nav-link", id: "navlig-lane3-tab", "data-toggle": "tab", href: "#navlig-lane3", role: "tab", "aria-controls": "navlig-lane3", "aria-selected": "false" },
          "Lane 3"
        ),
        React.createElement(
          "a",
          { className: "nav-item nav-link", id: "navlig-lane4-tab", "data-toggle": "tab", href: "#navlig-lane4", role: "tab", "aria-controls": "navlig-lane4", "aria-selected": "false" },
          "Lane 4"
        )
      )
    ),
    React.createElement(
      "div",
      { className: "tab-content", id: "nav-ligContent" },
      props.lig_tabs
    )
  );
}

function SummaryTable(props) {
  return React.createElement(
    "div",
    { className: "tab-pane fade show active", id: "summary", role: "tabpanel", "aria-labelledby": "summary-tab" },
    React.createElement(
      "div",
      { className: "d-flex justify-content-between flex-wrap flex-md-nowrap align-items-center pt-3 pb-2 mb-3" },
      React.createElement(
        "h1",
        { className: "h3", id: "summary-name" },
        "Summary Table"
      )
    ),
    React.createElement(
      "table",
      { className: "table table-hover" },
      React.createElement(
        "thead",
        null,
        React.createElement(
          "tr",
          null,
          React.createElement("th", { scope: "col" }),
          props.lane_stats.map(function (item, index) {
            return React.createElement(TitleRow, { key: index, samp: item.Lane });
          })
        )
      ),
      React.createElement(
        "tbody",
        null,
        React.createElement(
          "tr",
          null,
          React.createElement(
            "th",
            { scope: "row" },
            "Total input reads"
          ),
          props.lane_stats.map(function (item, index) {
            return React.createElement(RegRow, { key: index, val: item.tot_inp_reads });
          })
        ),
        React.createElement(
          "tr",
          null,
          React.createElement(
            "th",
            { scope: "row" },
            "Total passed reads"
          ),
          props.lane_stats.map(function (item, index) {
            return React.createElement(RegRow, { key: index, val: item.tot_pass_reads });
          })
        ),
        React.createElement(
          "tr",
          null,
          React.createElement(
            "th",
            { scope: "row" },
            "Pass percentage"
          ),
          props.lane_stats.map(function (item, index) {
            return React.createElement(RegRow, { key: index, val: item.pass_perc });
          })
        ),
        React.createElement(
          "tr",
          null,
          React.createElement(
            "th",
            { scope: "row", title: "Percent of reads where one of the barcode pieces was uncorrectable" },
            "Percent bad barcodes"
          ),
          props.lane_stats.map(function (item, index) {
            return React.createElement(RegRow, { key: index, val: item.perc_uncorr });
          })
        ),
        React.createElement(
          "tr",
          null,
          React.createElement(
            "th",
            { scope: "row", title: "Percent of reads where the RT barcode was not specified in the sample sheet" },
            "Percent invalid RT"
          ),
          props.lane_stats.map(function (item, index) {
            return React.createElement(RegRow, { key: index, val: item.perc_inval_rt });
          })
        ),
        React.createElement(
          "tr",
          null,
          React.createElement(
            "th",
            { scope: "row", title: "Percent of reads where the PCR barcodes were not one of the matched pairs provided" },
            "Percent PCR mismatch"
          ),
          props.lane_stats.map(function (item, index) {
            return React.createElement(RegRow, { key: index, val: item.perc_pcr_mismatch });
          })
        )
      )
    )
  );
}

function Header(props) {
  return React.createElement(
    "nav",
    { className: "navbar navbar-expand-md sticky-top navbar-light", style: { backgroundColor: "#e3f2fd" } },
    React.createElement(
      "div",
      { className: "navbar-collapse collapse w-100 order-1 order-md-0 dual-collapse2" },
      React.createElement(
        "ul",
        { className: "navbar-nav mr-auto" },
        React.createElement("img", { src: "img/bbi_icon.png", height: "70", className: "d-inline-block align-top", alt: "" })
      )
    ),
    React.createElement(
      "div",
      { className: "mx-auto order-0" },
      React.createElement(
        "a",
        { className: "navbar-brand mx-auto", href: "#" },
        "Demultiplexing ",
        props.run_name,
        " QC Dashboard"
      )
    ),
    React.createElement("div", { className: "navbar-collapse collapse w-100 order-3 dual-collapse2" })
  );
}

function DemuxPage(props) {
  return React.createElement(
    "span",
    null,
    React.createElement(Header, { run_name: props.run_name }),
    React.createElement(
      "div",
      { className: "container-fluid" },
      React.createElement(
        "div",
        { className: "row" },
        React.createElement(
          "nav",
          { className: "col-md-2 d-none d-md-block bg-light sidebar" },
          React.createElement(
            "div",
            { className: "sidebar-sticky" },
            React.createElement(
              "div",
              { className: "nav flex-column nav-pills", id: "v-pills-tab", role: "tablist", "aria-orientation": "vertical" },
              React.createElement(
                "a",
                { className: "nav-link active", id: "summary-tab", "data-toggle": "pill", href: "#summary", role: "tab", "aria-controls": "summary", "aria-selected": "true" },
                "Summary Table"
              ),
              React.createElement(
                "a",
                { className: "nav-link", id: "rt-tab", "data-toggle": "pill", href: "#rt", role: "tab", "aria-controls": "rt", "aria-selected": "false" },
                "RT Barcodes"
              ),
              React.createElement(
                "a",
                { className: "nav-link", id: "pcr-tab", "data-toggle": "pill", href: "#pcr", role: "tab", "aria-controls": "pcr", "aria-selected": "false" },
                "PCR Barcodes"
              ),
              props.level == 3 ? React.createElement(
                "a",
                { className: "nav-link", id: "lig-tab", "data-toggle": "pill", href: "#lig", role: "tab", "aria-controls": "lig", "aria-selected": "false" },
                "Ligation Barcodes"
              ) : ''
            )
          )
        ),
        React.createElement(
          "main",
          { role: "main", className: "col-md-9 ml-sm-auto col-lg-10 px-4", style: { paddingTop: "15px" } },
          React.createElement(
            "div",
            { className: "tab-content", id: "nav-tabContent" },
            React.createElement(SummaryTable, { lane_stats: props.lane_stats }),
            React.createElement(RTBarcodes, { rt_tabs: props.rt_tabs, bad_wells: props.bad_wells }),
            React.createElement(PCRBarcodes, { pcr_tabs: props.pcr_tabs }),
            props.level == 3 ? React.createElement(LigBarcodes, { lig_tabs: props.lig_tabs }) : ""
          )
        )
      )
    )
  );
}

/********** Generate pages **********/

var lane_list = run_data['lane_list'];
var lane_stats = run_data['lane_stats'];
var run_name = run_data['run_name'];
var plate_list = run_data['plate_list'];
var pcr_combo_list = run_data['pcr_combo_list'];
var lig_combo_list = run_data['lig_combo_list'];
var level = run_data['level'];
var bad_wells = run_data['bad_wells'];
var bad_wells_barcodes = run_data['bad_wells_barcodes'];
var include_norm = run_data['include_norm'];

var LigTabs = lane_list.map(function (lane, index) {
  return React.createElement(LigTab, { key: index, className: "tab-pane fade", lane: lane });
});

var RTTabs = lane_list.map(function (lane, index) {
  return React.createElement(RTTab, { key: index, className: "tab-pane fade", lane: lane, bad_wells: bad_wells, bad_wells_barcodes: bad_wells_barcodes, norm: include_norm });
});

var PCRLaneTabs = lane_list.map(function (lane, index) {
  return React.createElement(PCRTab, { key: index, className: "tab-pane fade", lane: lane });
});

ReactDOM.render(React.createElement(DemuxPage, { rt_tabs: RTTabs, pcr_tabs: PCRLaneTabs, lig_tabs: LigTabs, lane_stats: lane_stats, level: level, run_name: run_name, bad_wells: bad_wells }), document.getElementById('demux_page'));