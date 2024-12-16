// function handle_submit() {

//     var request = new XMLHttpRequest();
//     // console.log(button_european.value);
//     // console.log("https://localhost:7007/" + button_european.value);
//     request.open("G", "http://localhost:5115/Trade");
//     request.setRequestHeader("Content-Type", "application/json");

//     request.onreadystatechange = function () {
//         if (request.readyState === 4 && request.status === 200) {
//             var output = JSON.parse(request.responseText);
//             console.log(output);
//         }
//     };

//     var form_data = new FormData(document.forms.form_european);

//     var object = {};
//     form_data.forEach(function(value, key) {
//         if (value === 'True') {
//             value = true;
//         } 

//         if (value === 'False') {
//             value = false;
//         }

//         object[key] = value;
//     });
//     var json = JSON.stringify(object);

//     console.log(json)
//     request.send(json);

// }

// const button_european = document.getElementById("valuation_button");

// button_european.addEventListener("click", handle_submit);

// var radios = document.forms["option_select"].elements["option_name"];
// for(var i = 0, max = radios.length; i < max; i++) {
//     radios[i].onclick = function() {
//         var option_display = document.getElementById("option_select_display");
//         button_european.value = this.value;
//         option_display.innerHTML = button_european.value + " Option";
//         button_european.innerHTML = "Simulate " + this.value;
//         // myFunction();
//         console.log(button_european.value);

//         var IsCall_button = document.getElementById("IsCall_input");
//         if (button_european.value == "Range" & IsCall_button.style.display == "block") {
//             IsCall_button.style.display = "none";
//         }
//         else {
//             IsCall_button.style.display = "block";
//         }

//         var strike_button = document.getElementById("strike_input");
//         if (button_european.value == "Range" & strike_button.style.display == "block") {
//             strike_button.style.display = "none";
//         }
//         else {
//             strike_button.style.display = "block";
//         }

//         var payout_button = document.getElementById("payout_input");
//         if (button_european.value == "Digital" & payout_button.style.display == "none") {
//             payout_button.style.display = "block";
//         }
//         else {
//             payout_button.style.display = "none";
//         }

//         var barrierlevel_button = document.getElementById("barrierlevel_input");
//         if (button_european.value == "Barrier" & barrierlevel_button.style.display == "none") {
//             barrierlevel_button.style.display = "block";
//         }
//         else {
//             barrierlevel_button.style.display = "none";
//         }

//         var knocktype_button = document.getElementById("knocktype_input");
//         if (button_european.value == "Barrier" & knocktype_button.style.display == "none") {
//             knocktype_button.style.display = "block";
//         }
//         else {
//             knocktype_button.style.display = "none";
//         }
//     }
// }

var all_european = [];
var all_digital = [];
var all_asian = [];
var all_range = [];
var all_lookback = [];
var all_barrier = [];
let all_results = [];

var one_day = 1000 * 60 * 60 * 24;
var today = new Date();

function handle_get() {

    var request = new XMLHttpRequest();
    var trade = new XMLHttpRequest();

    // result to empty list before running requests
    all_european = [];
    all_digital = [];
    all_asian = [];
    all_range = [];
    all_lookback = [];
    all_barrier = [];
    all_results = [];

    request.open("Get", "http://localhost:5115/OptionsData");
    request.setRequestHeader("Content-Type", "application/json");

    trade.open("Get", "http://localhost:5115/Trade");
    trade.setRequestHeader("Content-Type", "application/json");

    trade.onreadystatechange = function () {
        if (trade.readyState === 4 && trade.status === 200) {
            var tradeoutput = JSON.parse(trade.responseText);
            console.log("trade output", tradeoutput);

            request.onreadystatechange = function () {
                if (request.readyState === 4 && request.status === 200) {
                    var output = JSON.parse(request.responseText);
                    console.log("JS output", output);
        
                    output.forEach(x => {
                        var mydict = {};
                        var sim_params = handle_simulation_params();
                        var expiry = new Date(x[Object.keys(x)[0]].expiration_date);
                        
                        var tenor = Math.round(expiry.getTime() - today.getTime()) / (one_day) / 365; 

                        tradeoutput.forEach(y => {
                            if (y["financialinstrumentid"] == x[Object.keys(x)[0]]["id"])
                            {
                                sim_params["S"] = y["trade_price"];
                            }
                        });
                        sim_params["Strike"] = x[Object.keys(x)[0]]["strike"];
                        sim_params["T"] = tenor;
                        sim_params["id"] = x[Object.keys(x)[0]]["id"];
        
                        if (Object.keys(x)[0] == "European") {
                            all_european.push(sim_params);
                        }
                        else if (Object.keys(x)[0] == "Digital") {
                            sim_params["Payout"] =  x[Object.keys(x)[0]]["payout"];
                            all_digital.push(sim_params);
                        }
                        else if (Object.keys(x)[0] == "Asian") {
                            all_asian.push(sim_params);
                        }
                        else if (Object.keys(x)[0] == "Lookback") {
                            all_lookback.push(sim_params);
                        }
                        else if (Object.keys(x)[0] == "Range") {
                            all_range.push(sim_params);
                        }
                        else if (Object.keys(x)[0] == "Barrier") {
                            sim_params["barrier_level"] =  x[Object.keys(x)[0]]["barrier_Level"];
                            sim_params["knock_type"] =  x[Object.keys(x)[0]]["knock_Type"];
                            all_barrier.push(sim_params);
                        }
                    });
        
                    console.log(all_european);
                    console.log(all_digital);
                    console.log(all_asian);
                    console.log(all_range);
                    console.log(all_lookback);
                    console.log(all_barrier);
        
                    all_european.forEach(x => {
                        var request = new XMLHttpRequest();
                        request.open("POST", "http://localhost:5115/European");
                        request.setRequestHeader("Content-Type", "application/json");

                        request.onreadystatechange = function () {
                            if (request.readyState === 4 && request.status === 200) {
                                var output = JSON.parse(request.responseText);
                                tradeoutput.forEach(y => {
                                    if (y["financialinstrumentid"] == x["id"])
                                    {
                                        output["quantity"] = y["quantity"];
                                        output["financialinstrumentid"] = y["financialinstrumentid"];
                                        output["trade_price"] = y["trade_price"];
                                        output["trade_id"] = y["id"];
                                    }
                                });
                                console.log("euro all res", all_results);
                                all_results.push(output);
                            }
                        };
        
                        var json = JSON.stringify(x);
        
                        console.log("json", json)
                        request.send(json);
                    });
        
                    all_asian.forEach(x => {
                        var request = new XMLHttpRequest();
                        request.open("POST", "http://localhost:5115/Asian");
                        request.setRequestHeader("Content-Type", "application/json");
        
                        request.onreadystatechange = function () {
                            if (request.readyState === 4 && request.status === 200) {
                                var output = JSON.parse(request.responseText);
                                tradeoutput.forEach(y => {
                                    if (y["financialinstrumentid"] == x["id"])
                                    {
                                        output["quantity"] = y["quantity"];
                                        output["financialinstrumentid"] = y["financialinstrumentid"];
                                        output["trade_price"] = y["trade_price"];
                                        output["trade_id"] = y["id"];
                                    }
                                });
                                all_results.push(output);
                            }
                        };
        
                        var json = JSON.stringify(x);
        
                        console.log("json", json)
                        request.send(json);
                    });
        
                    all_digital.forEach(x => {
                        var request = new XMLHttpRequest();
                        request.open("POST", "http://localhost:5115/Digital");
                        request.setRequestHeader("Content-Type", "application/json");
        
                        request.onreadystatechange = function () {
                            if (request.readyState === 4 && request.status === 200) {
                                var output = JSON.parse(request.responseText);
                                tradeoutput.forEach(y => {
                                    if (y["financialinstrumentid"] == x["id"])
                                    {
                                        output["quantity"] = y["quantity"];
                                        output["financialinstrumentid"] = y["financialinstrumentid"];
                                        output["trade_price"] = y["trade_price"];
                                        output["trade_id"] = y["id"];
                                    }
                                });
                                all_results.push(output);
                            }
                        };
        
                        var json = JSON.stringify(x);
        
                        console.log("json", json)
                        request.send(json);
                    });
        
                    all_range.forEach(x => {
                        var request = new XMLHttpRequest();
                        request.open("POST", "http://localhost:5115/Range");
                        request.setRequestHeader("Content-Type", "application/json");
        
                        request.onreadystatechange = function () {
                            if (request.readyState === 4 && request.status === 200) {
                                var output = JSON.parse(request.responseText);
                                tradeoutput.forEach(y => {
                                    if (y["financialinstrumentid"] == x["id"])
                                    {
                                        output["quantity"] = y["quantity"];
                                        output["financialinstrumentid"] = y["financialinstrumentid"];
                                        output["trade_price"] = y["trade_price"];
                                        output["trade_id"] = y["id"];
                                    }
                                });
                                all_results.push(output);
                            }
                        };
        
                        var json = JSON.stringify(x);
        
                        console.log("json", json)
                        request.send(json);
                    });
        
                    all_lookback.forEach(x => {
                        var request = new XMLHttpRequest();
                        request.open("POST", "http://localhost:5115/Lookback");
                        request.setRequestHeader("Content-Type", "application/json");
        
                        request.onreadystatechange = function () {
                            if (request.readyState === 4 && request.status === 200) {
                                var output = JSON.parse(request.responseText);
                                tradeoutput.forEach(y => {
                                    if (y["financialinstrumentid"] == x["id"])
                                    {
                                        output["quantity"] = y["quantity"];
                                        output["financialinstrumentid"] = y["financialinstrumentid"];
                                        output["trade_price"] = y["trade_price"];
                                        output["trade_id"] = y["id"];
                                    }
                                });
                                all_results.push(output);
                            }
                        };
        
                        var json = JSON.stringify(x);
        
                        console.log("json", json)
                        request.send(json);
                    });
        
                    all_barrier.forEach(x => {
                        var request = new XMLHttpRequest();
                        request.open("POST", "http://localhost:5115/barrier");
                        request.setRequestHeader("Content-Type", "application/json");
        
                        request.onreadystatechange = function () {
                            if (request.readyState === 4 && request.status === 200) {
                                var output = JSON.parse(request.responseText);
                                tradeoutput.forEach(y => {
                                    if (y["financialinstrumentid"] == x["id"])
                                    {
                                        output["quantity"] = y["quantity"];
                                        output["financialinstrumentid"] = y["financialinstrumentid"];
                                        output["trade_price"] = y["trade_price"];
                                        output["trade_id"] = y["id"];
                                    }
                                });
                                all_results.push(output);
                            }
                        };
        
                        var json = JSON.stringify(x);
        
                        console.log("json", json)
                        request.send(json);
                    });


                    // combine all info into a dict matching
                    console.log("allres", all_results);
    
                }
            }
        }
    }

    trade.send();
    request.send();
}

function update_table() {

    console.log("all_results update", all_results.length);

    var mytable = document.getElementById("tradestable");
    let rowcount = mytable.rows.length;

    for (let i = 1; i < rowcount; i++) {
        mytable.deleteRow(-1);
    }

    for (let i = 0; i < all_results.length; i++) {
        var row = mytable.insertRow(mytable.rows.length);
        
        // trade id
        var mytradeid = row.insertCell();
        mytradeid.innerHTML = all_results[i].trade_id;

        var myfinancialinstrumentid = row.insertCell();
        myfinancialinstrumentid.innerHTML = all_results[i].financialinstrumentid;

        var myquantity = row.insertCell();
        myquantity.innerHTML = all_results[i].quantity;

        var mytradeprice = row.insertCell();
        mytradeprice.innerHTML = all_results[i].trade_price.toFixed(2);

        var mysimprice = row.insertCell();
        mysimprice.innerHTML = all_results[i].price.toFixed(4);

        var myse = row.insertCell();
        myse.innerHTML = all_results[i].se.toFixed(4);

        var mydelta = row.insertCell();
        mydelta.innerHTML = all_results[i].delta.toFixed(2);

        var mygamma = row.insertCell();
        mygamma.innerHTML = all_results[i].gamma.toFixed(4);

        var myvega = row.insertCell();
        myvega.innerHTML = all_results[i].vega.toFixed(4);

        var mytheta = row.insertCell();
        mytheta.innerHTML = all_results[i].theta.toFixed(4);

        var myrho = row.insertCell();
        myrho.innerHTML = all_results[i].rho.toFixed(4);

        var mypnl = row.insertCell();
        mypnl.innerHTML = (all_results[i].quantity * (all_results[i].price - all_results[i].trade_price)).toFixed(2);

    }
}

const button_get = document.getElementById("valuation_refresh");
button_get.addEventListener("click", update_table);

function handle_simulation_params() {
    var form_data = new FormData(document.forms.form_params);

    var object = {};
    form_data.forEach(function(value, key) {
        if (value === 'True') {
            value = true;
        } 

        if (value === 'False') {
            value = false;
        }

        object[key] = value;
    });
    object["r"] = Number(object["r"]);
    object["div"] = Number(object["div"]);
    object["sig"] = Number(object["sig"]);
    object["N"] = Number(object["N"]);
    object["M"] = Number(object["M"]);


    return object;
}

const button = document.getElementById("valuation_button");
button.addEventListener("click", handle_get);