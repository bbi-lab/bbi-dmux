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
        props.rt_tab_head
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
    { className: props.className, id: "navpcr-lane" + props.lane, role: "tabpanel",
      "aria-labelledby": "navpcr-lane" + props.lane + "-tab" },
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
    React.createElement("img", { src: "img/L00" + props.lane + "_" + props.plate_name + ".pcr_plate.png", width: "50%",
      className: "rounded mx-auto d-block", alt: "..." })
  );
}

function PCRBarcodes(props) {
  var lane_list = props.lane_stats.map(function (item) {
    return item.Lane;
  });
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
    Object.keys(props.pcr_combo_list).length === 0 ? React.createElement(
      "span",
      null,
      React.createElement(
        "h4",
        null,
        "Well combination read counts:"
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
            React.createElement(
              "th",
              { scope: "col" },
              "p5 well"
            ),
            React.createElement(
              "th",
              { scope: "col" },
              "p7 well"
            ),
            lane_list.map(function (item, index) {
              return React.createElement(TitleRow, { key: index, samp: item });
            })
          )
        ),
        React.createElement(
          "tbody",
          null,
          props.pcr_well_info.map(function (item, index) {
            return React.createElement(
              "tr",
              { key: index },
              React.createElement(
                "th",
                { scope: "row" },
                item.p5
              ),
              React.createElement(
                "th",
                { scope: "row" },
                item.p7
              ),
              lane_list.map(function (lane, index) {
                return React.createElement(RegRow, { key: index, val: item[lane] });
              })
            );
          })
        )
      )
    ) : React.createElement(
      "span",
      null,
      React.createElement(
        "nav",
        null,
        React.createElement(
          "div",
          { className: "nav nav-tabs", id: "navpcr-tab", role: "tablist" },
          props.pcr_tab_head
        )
      ),
      React.createElement(
        "div",
        { className: "tab-content", id: "nav-pcrContent" },
        props.pcr_tabs
      )
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

function CodeChunk(props) {
  return React.createElement(
    "pre",
    { style: { paddingLeft: '20px' } },
    React.createElement(
      "code",
      null,
      '\n' + props.text + '\n\n'
    )
  );
}

function LogTab(props) {
  return React.createElement(
    "div",
    { className: props.className, id: "navlog-lane" + props.lane, role: "tabpanel", "aria-labelledby": "navlog-lane" + props.lane + "-tab" },
    React.createElement(CodeChunk, { text: props.log })
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
        props.lig_tab_head
      )
    ),
    React.createElement(
      "div",
      { className: "tab-content", id: "nav-ligContent" },
      props.lig_tabs
    )
  );
}

function RecoveryLog(props) {
  return React.createElement(
    "div",
    { className: "tab-pane fade", id: "log", role: "tabpanel", "aria-labelledby": "log-tab" },
    React.createElement(
      "div",
      { className: "d-flex justify-content-between flex-wrap flex-md-nowrap align-items-center pt-3 pb-2 mb-3 border-bottom" },
      React.createElement(
        "h1",
        { className: "h3", id: "log-name" },
        "Recovery Summary"
      )
    ),
    React.createElement(
      "nav",
      null,
      React.createElement(
        "div",
        { className: "nav nav-tabs", id: "navlig-tab", role: "tablist" },
        props.log_tab_head
      )
    ),
    React.createElement(
      "div",
      { className: "tab-content", id: "nav-ligContent" },
      props.log_tabs
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
              ) : '',
              React.createElement(
                "a",
                { className: "nav-link", id: "log-tab", "data-toggle": "pill", href: "#log", role: "tab", "aria-controls": "log", "aria-selected": "false" },
                "Recovery Summary"
              )
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
            React.createElement(RTBarcodes, { rt_tabs: props.rt_tabs, bad_wells: props.bad_wells, rt_tab_head: props.rt_tab_head }),
            React.createElement(PCRBarcodes, { pcr_tabs: props.pcr_tabs, lane_stats: props.lane_stats,
              pcr_combo_list: props.pcr_combo_list, pcr_well_info: props.pcr_well_info,
              pcr_tab_head: props.pcr_tab_head }),
            props.level == 3 ? React.createElement(LigBarcodes, { lig_tabs: props.lig_tabs, lig_tab_head: props.lig_tab_head }) : "",
            React.createElement(RecoveryLog, { log_tabs: props.log_tabs, log_tab_head: props.log_tab_head })
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
var pcr_well_info = run_data['pcr_well_info'];

var LigTabs = lane_list.map(function (lane, index) {
  return React.createElement(LigTab, { key: index, className: "tab-pane fade", lane: lane });
});

var LogTabs = lane_list.map(function (lane, index) {
  return React.createElement(LogTab, { key: index, className: "tab-pane fade", lane: lane, log: log_data[lane] });
});

var RTTabs = lane_list.map(function (lane, index) {
  return React.createElement(RTTab, { key: index, className: "tab-pane fade", lane: lane, bad_wells: bad_wells, bad_wells_barcodes: bad_wells_barcodes, norm: include_norm });
});

var PCRLaneTabs = lane_list.map(function (lane, index) {
  return Object.keys(pcr_combo_list).length === 0 ? "" : React.createElement(PCRTab, { key: index, className: "tab-pane fade", lane: lane });
});

var LigTabHead = lane_list.map(function (lane, index) {
  return React.createElement(
    "a",
    { key: index, className: "nav-item nav-link", id: "navlig-lane" + lane + "-tab", "data-toggle": "tab", href: "#navlig-lane" + lane, role: "tab", "aria-controls": "navlig-lane" + lane, "aria-selected": "false" },
    "Lane " + lane
  );
});

var RTTabHead = lane_list.map(function (lane, index) {
  return React.createElement(
    "a",
    { key: index, className: "nav-item nav-link", id: "navrt-lane" + lane + "-tab", "data-toggle": "tab", href: "#navrt-lane" + lane, role: "tab", "aria-controls": "navrt-lane" + lane, "aria-selected": "false" },
    "Lane " + lane
  );
});

var PCRTabHead = lane_list.map(function (lane, index) {
  return React.createElement(
    "a",
    { key: index, className: "nav-item nav-link", id: "navpcr-lane" + lane + "-tab", "data-toggle": "tab", href: "#navpcr-lane" + lane, role: "tab", "aria-controls": "navpcr-lane" + lane, "aria-selected": "false" },
    "Lane " + lane
  );
});

var LogTabHead = lane_list.map(function (lane, index) {
  return React.createElement(
    "a",
    { key: index, className: "nav-item nav-link", id: "navlog-lane" + lane + "-tab", "data-toggle": "tab", href: "#navlog-lane" + lane, role: "tab", "aria-controls": "navlog-lane" + lane, "aria-selected": "false" },
    "Lane " + lane
  );
});

ReactDOM.render(React.createElement(DemuxPage, { rt_tabs: RTTabs, pcr_tabs: PCRLaneTabs, pcr_well_info: pcr_well_info,
  pcr_combo_list: pcr_combo_list, lig_tabs: LigTabs, lane_stats: lane_stats,
  level: level, run_name: run_name, bad_wells: bad_wells, rt_tab_head: RTTabHead,
  lig_tab_head: LigTabHead, pcr_tab_head: PCRTabHead, log_tabs: LogTabs, log_tab_head: LogTabHead }), document.getElementById('demux_page'));