const button = document.getElementById("derivatives_submit");
document.onload = handle_get();

function handle_submit() {

    var request = new XMLHttpRequest();
    // console.log(button_european.value);
    // console.log("https://localhost:7007/" + button_european.value);
    request.open("POST", "http://localhost:5115/" + button.value +"PostData");
    request.setRequestHeader("Content-Type", "application/json");

    request.onreadystatechange = function () {
        if (request.readyState === 4 && request.status === 200) {
            var output = JSON.parse(request.responseText);
            console.log(output);
        }
    };

    var form_data = new FormData(document.forms.form_derivatives);

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
    var json = JSON.stringify(object);

    console.log(json)
    request.send(json);
}


button.addEventListener("click", handle_submit);

var radios = document.forms["option_select"].elements["option_name"];
for(var i = 0, max = radios.length; i < max; i++) {
    radios[i].onclick = function() {
        var option_display = document.getElementById("option_select_display");
        button.value = this.value;
        option_display.innerHTML = button.value + " Option";
        button.innerHTML = "Submit " + this.value;
        button_get.innerHTML = "Refresh " + this.value;
        // myFunction();
        console.log(button.value);

        var IsCall_button = document.getElementById("IsCall_input");
        if (button.value == "Range" & IsCall_button.style.display == "block") {
            IsCall_button.style.display = "none";
        }
        else {
            IsCall_button.style.display = "block";
        }

        var strike_button = document.getElementById("strike_input");
        if (button.value == "Range" & strike_button.style.display == "block") {
            strike_button.style.display = "none";
            strike_button.disabled = true;
            strike_button.value = null;
        }
        else {
            strike_button.style.display = "block";
        }

        var payout_button = document.getElementById("payout_input");
        if (button.value == "Digital" & payout_button.style.display == "none") {
            payout_button.style.display = "block";
        }
        else {
            payout_button.style.display = "none";
        }

        var barrierlevel_button = document.getElementById("barrierlevel_input");
        if (button.value == "Barrier" & barrierlevel_button.style.display == "none") {
            barrierlevel_button.style.display = "block";
        }
        else {
            barrierlevel_button.style.display = "none";
        }

        var knocktype_button = document.getElementById("knocktype_input");
        if (button.value == "Barrier" & knocktype_button.style.display == "none") {
            knocktype_button.style.display = "block";
        }
        else {
            knocktype_button.style.display = "none";
        }

        handle_get();
    }
}

function handle_get() {

    var request = new XMLHttpRequest();
    request.open("Get", "http://localhost:5115/"+ button.value + "GetData");
    request.setRequestHeader("Content-Type", "application/json");

    request.onreadystatechange = function () {
        if (request.readyState === 4 && request.status === 200) {
            var output = JSON.parse(request.responseText);
            console.log(output);
            var mytable = document.getElementById("derivativestable");
            let rowcount = mytable.rows.length;
            for (let i = 1; i < rowcount; i++) {
                mytable.deleteRow(-1);
            }
            for (let i = 0; i < output.length; i++) {
                var row = mytable.insertRow(mytable.rows.length);

                var idCol = row.insertCell();
                idCol.innerHTML = output[i].id;

                var mymarketid = row.insertCell();
                mymarketid.innerHTML = output[i].marketid;

                var myunderlyingid = row.insertCell();
                myunderlyingid.innerHTML = output[i].underlyingid;

                var mydate = row.insertCell();
                mydate.innerHTML = output[i].expiration_date;

                var mystrike = row.insertCell();
                mystrike.innerHTML = output[i].strike;

                var mycall = row.insertCell();
                mycall.innerHTML = output[i].is_Call;

                var mypayout = row.insertCell();
                mypayout.innerHTML = output[i].payout;

                var mybarrierlevel = row.insertCell();
                mybarrierlevel.innerHTML = output[i].barrier_Level;

                var myknock = row.insertCell();
                myknock.innerHTML = output[i].knock_Type;
            }
        }
    };
    console.log("JS Get");
    request.send();
}

const button_get = document.getElementById("derivatives_refresh");
button_get.addEventListener("click", handle_get);